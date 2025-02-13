% -------------------------------------------------------
%
%    preprocess12 - function to execute initial preprocessing steps
%
%    Created:         Daniela Zoeller      (30.10.2015)
%    Based on: preprocess.m (Jonas Richiardi and Giulia Preti)
%    Last modified:   Daniela Zoeller      (15.03.2016)
%
% ------------------------------------------------------
%
% Modification preprocess.m:
%   30.10.2015: - run newSegment with Dartel option
%               - no deformation task and registration to MNI
%   15.03.2016: - modification of functions and jobs for SPM12

function preprocess12(functPath, structPath,prefix_struct, varargin)
% Preprocessing pipeline for fMRI data (bias field correction, functional
% realignment, functional-structural coregistration, segmentation,
% normalisation, atlasing)
%
% IN:
%   functPath: full path to raw fatlasFileunctional data
%   structPath: full path to raw structural data
%   varargin:
%       1)  the processing chain to use, e.g.
%       'procChain',{'realign','QC','coregister','newSegment','label'}
%       2) other arguments for quality control
%       'QCcoef', 1.5
%       3) high-pass filtering
%       'highpass', 0
%
% REQUIREMENTS
% - SPM8
% - IBASPM toolbox
% 
% INSTALLATION
% - SPM8 must be on your matlab path
%
% ACKNOWLEDGEMENTS
% This file borrows liberally from the following sources
% - Rik Henson's example script for SPM 5 (MRC Cognition and Brain Sciences Unit)
% http://en.wikibooks.org/wiki/SPM-Example_batch_script
% - Various functions in the IBASPM toolbox, from Yasser Aleman-Gomez and
% colleages at the Cuba Neurosciences Center
%
% v1.0 Jonas Richiardi / Medical Image Processing Laboratory
% - initial release for SPM5
% v2.0 Manuel W??thrich
% - SPM8 support
% - use spm jobmanager
% - support for newSegment and unifiedSegmentation
% - port IBASPM auto_labelling to SPM8
% - basic QC capability
% v2.0.1 Jonas Richiardi
% - cross-platformish
% - path settings fixed
% v2.0.1b Jonas Richiardi
% - removed ~ ignore code and changed to 4-output forme for filepart for 
% backwards compatibility with pre-7.9 versions e.g. as in Greedy SCC cluster
% - corrected path setting code
% v2.0.2 June 2011 Jonas Richiardi
% - supports smoothing
% - supports coregistration direction setting (FtoS or StoF)
% v2.0.3 June 2011 JR
% - supports custom atlases
% v2.1 Oct 2011 JR
% - fixes a critical bug introduced in v2.0 - wrong interpolation
% parameters were used in new_Auto_Labelling.m
% v2.1.1 Dec 2011 JR
% - stricter dir-and-file checking in pp_loadVolumes
% v2.1.2 Dec 2011 JR
% - add version tagging to tasksDone structure to help traceability
% - bugfix and better existence type checking in pp_loadVolume
% - support for Matlab 7.13 (fileparts syntax)

DEBUGMODE=true; % set this to true to help diagnosing problems

RestingState_mode = 1;

%% check and set necessary paths
if exist('spm.m','file')~=2
    error('SPM seems not to be on your matlab path.');
else
    spmLoc=spm('Dir');
    myPath=path();
    % add 'config' SPM dir to the path
    spmConfigLoc=[spmLoc filesep 'config'];
    if exist(spmConfigLoc,'dir')~=7
        error(['SPM installation seems to be missing a config dir at '...
            spmConfigLoc ', quitting.']);
    else
        addpath(spmConfigLoc);
    end
end

if (DEBUGMODE==true)
    path
end

%% initialize and load data
disp('** STARTING PROCESSING');
Tstart = clock;
spm('defaults', 'fMRI');
pp12_loadVolumes;
% generate full path of this file
thisLocation=which('preprocess12.m');
% max-compatibility version of fileparts (old releases don't have ~)
if verLessThan('matlab', '7.13.0')
    [jobsParentPath, quux1, quux2, quux3] = fileparts(thisLocation);
else
    [jobsParentPath, quux1, quux2] = fileparts(thisLocation);
end


%% reset of origin
if tasksTodo.reset
    % -------------------- reset of functional images -------------------------------------
    
    cd(functPath);  
    P = dir('*_*');
    %[P,sts] = spm_select('List',functPath,'dir','_*');
    %if ~sts, return; else P = cellstr(P); end
    spm_progress_bar('Init',numel(P),'Resetting orientations',...
        'Images Complete');
    for i=1:numel(P)
        %V    = spm_vol(P{i});
        V    = spm_vol(P(i).name);
        M    = V.mat;
        vox  = sqrt(sum(M(1:3,1:3).^2));
        if det(M(1:3,1:3))<0, vox(1) = -vox(1); end
        orig = (V.dim(1:3)+1)/2;
        off  = -vox.*orig;
        M    = [vox(1) 0      0      off(1)
            0      vox(2) 0      off(2)
            0      0      vox(3) off(3)
            0      0      0      1];
        spm_get_space(P(i).name,M);
        spm_progress_bar('Set',i);
    end
    spm_progress_bar('Clear');
end

%---------------- update tasksTodo -------------------------
    tasksDone.reset = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
    
    
%% realign

if  tasksTodo.realign
    % -------------------- realign -------------------------------------
    cd(functPath);
%     if exist(alignFolder,'dir'); rmdir(alignFolder,'s'); end;
%     mkdir(alignFolder);
    jobs = {fullfile(jobsParentPath,'jobs','align_job.mat')};
    spm_jobman('initcfg');
    spm_jobman('serial', jobs, '', functFiles);

    % ----------------- move files ----------------------
%     system(strcat('find -name "rf*" -exec mv {} ',functPath,' \;'));
%     movelist = dir(fullfile(functPath,'rf*'));  
%     ind = [1, round(length(movelist)/3), round(2*length(movelist)/3), length(movelist)];
%     for j = 1:3
%         movelist2 = movelist(ind(j):ind(j+1));
%         folders = cell(size(movelist2));
%         for i = 1:size(folders,1)
%             folders{i} = fullfile(functPath,movelist2(i).name);
%         end
%         movefile(cellstr(folders),alignFolder);
%     end
    movelist = dir(fullfile(functPath,'r*_*.nii'));  
    for j = 1:length(movelist)
        movefile(fullfile(functPath, movelist(j).name),alignFolder);
    end
    
    movefile(fullfile(functPath,'rp*'),alignFolder);
    movefile(fullfile(functPath,'mean*'),alignFolder);
    copyfile(jobs{:}, jobFolder);
    
    %---------------- update tasksTodo -------------------------
    tasksDone.realign = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
    
    
end

%% QC
if tasksTodo.QC

    % -------------------- load ----------------------------------------
    cd(fullfile(alignFolder));
    file= dir('rp_*.txt');
    transData = load(fullfile(alignFolder, file.name));
    
    %functSize = size(spm_read_vols(spm_vol(functFiles{5})));
    functSize = 3;
    intens = zeros(size(functFiles,1),1);
    intensTop100 = zeros(size(functFiles,1),1);
    angle = zeros(size(transData,1),1);
    transl = zeros(size(transData,1),1);
    maxTransl = zeros(size(functFiles,1),1);
    corners = cell(8,1);
    
    for k = 0:7
        mask = dec2bin(k,3) == '1';
        corners{k+1} = (mask .* functSize)';
    end
    
    for i = 1:size(transData,1)
        transMatrix = spm_matrix(transData(i,:));
        quat = dcm2qua(transMatrix(1:3,1:3));
        angle(i) = 2*acos(quat(1))/(2*pi)*360;
        transl(i) = norm(transMatrix(1:3,4));
        for k = 1:8
            dist = norm(corners{k} - transMatrix(1:3,:) * [corners{k};1]);
            if dist > maxTransl(i), maxTransl(i) = dist; end
        end
        
        image = spm_read_vols(spm_vol(functFiles{i}));
        intens(i) = mean(image(:));
        array = sort(image(:),'descend');
        intensTop100(i) = mean(array(1:100));
    end    
   
    intens = (intens-mean(intens))/std(intens);
    intensTop100 = (intensTop100-mean(intensTop100))/std(intensTop100); 
    
    % ---------- compute outliers --------------------------------------
    x = (1:size(intens,1))';
    c = polyfit(x,intens,2);
    f = polyval(c,x);
    
    % ---------  on mean intensity -------------------------------------
    sortIntens = sort(intens);
    lQ = sortIntens(round(0.25*size(intens,1)));
    uQ = sortIntens(round(0.75*size(intens,1)));
    IQR = uQ - lQ;
    min = f + lQ - QCcoef*IQR;
    max = f + uQ + QCcoef*IQR;
    
    indexInt = [find(intens > max); find(intens < min)];
    
%     figure;
%     plot(intens,'b'); hold on;
%     plot(max,'g');
%     plot(min,'g');
%     plot(indexInt,intens(indexInt),'ro');hold off;
    
%    ---------- compute outliers --------------------------------------
%     x = (1:size(intensTop100,1))';
%     c = polyfit(x,intensTop100,2);
%     f = polyval(c,x);
    
    %----------- on mean intensity of top 100 voxels ------------------
%     sortIntensTop100 = sort(intensTop100);
%     lQ = sortIntensTop100(round(0.25*size(intensTop100,1)));
%     uQ = sortIntensTop100(round(0.75*size(intensTop100,1)));
%     IQR = uQ - lQ;
%     min = f + lQ - 1.5*IQR;
%     max = f + uQ + 1.5*IQR;
%     
%     indexTop100 = [find(intensTop100 > max); find(intensTop100 < min)];
%     
%     figure;
%     plot(intensTop100,'b'); hold on;
%     plot(max,'g');
%     plot(min,'g');
%     plot(indexTop100,intensTop100(indexTop100),'ro');hold off;
%     
    % ------------ apply ---------------------------------------------
    mask = zeros(size(intens,1),1);
    mask(indexInt) = 1;
    
    if highpass_bool
        mkdir(QCFolder);
        %alignFilenames = spm_select('List',alignFolder,['^r.*\.' volExt '$']);
        cd (alignFolder);
        files = dir(fullfile(alignFolder,'r*.nii'));
        alignFilename=zeros(length(files));
        for i=1:size(files)
            alignFilename(i)=files(i).name;
        end
        alignFiles=cell(size(alignFilenames,1),1);
        QCFiles=cell(size(alignFilenames,1),1);
        for f=1:size(alignFilenames,1)
            alignFiles{f}=fullfile(alignFolder,alignFilenames(f,:));
            QCFiles{f}=fullfile(QCFolder,alignFilenames(f,:));
        end

        for i = size(functFiles,1):-1:1
           if mask(i) == 0, copyfile(alignFiles{i},QCFolder);
           else QCFiles(i) = []; end
        end

        highpassFilter(char(QCFiles), 2, 0);

        for f=1:size(QCFiles,1)
            delete(QCFiles{f});
        end
    end
    
    
    
    % -------------- save artifacts ----------------------------------
    artifacts = struct('rotation', {angle}, 'translation', {transl},'maxTranslation',{maxTransl},...
        'intensity', {intens}, 'intensityTop100', {intensTop100},'mask',mask);
    save(fullfile(alignFolder, 'artifacts'), 'artifacts');
    
    %---------------- update tasksTodo -------------------------
    tasksDone.QC = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate
if  tasksTodo.coregister
    jobs = {fullfile(jobsParentPath,'jobs','coreg_job.mat')};
    inputs = cell(3, 1);
    structf = dir(fullfile(structPath,'rs*.nii'));
    structFile = fullfile(structPath,structf(1).name);
    coregDirection;
    if strcmp(coregDirection,'StoF')
            % coreg struct to func (change headers of struct file)
        %---------------- load meanfile -------------------------
        %functFilename = spm_select('List',alignFolder,['^mean.*\.*']);
        cd(alignFolder)
        files = dir(fullfile(alignFolder,'mean*'));
        meanFile=fullfile(alignFolder,files(1).name);
        %
            inputs{1} = cellstr(meanFile);
            inputs{2} =  cellstr(structFile); % Coreg: Estimate: Source Image - cfg_files
            display('\n Coreg StoF \n'); 
            inputs{3} = {''};
    else
        cd(alignFolder)
        files = dir(fullfile(alignFolder,'mean*'));
        meanFile=fullfile(alignFolder,files(1).name);
            % coreg funct to struct (change headers of func files)
            inputs{1} =  cellstr(structFile); % Coreg: Estimate: target Image
            inputs{2} = cellstr(meanFile);
            display('\n Coreg FtoS \n');
            file = dir('rp*.txt');
            tmp_others=cellstr(file.name);
            %tmp_others=cellstr(spm_select('List',alignFolder,'rp*.txt'));
            tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
            inputs{3} = tmp_others; % others (all other funcs)
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregister = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% coregister: estimate and and reslice
if  tasksTodo.coregReslice
    %structf = dir(fullfile(structPath,[prefix_struct,'*']));
    structf = dir(fullfile(structPath,'rs*.nii'));
    jobs = {fullfile(jobsParentPath,'jobs','coregReslice_job.mat')};
    inputs = cell(4, 1);
    
    functf = dir(fullfile(functPath,'realigned','r*.nii'));
    functf = struct2cell(functf); 
    functf = functf(1,:);
    functf_fullpath=cellfun(@(x) fullfile(functPath,'realigned',x), functf,'UniformOutput',false);

    
    if strcmp(coregDirection,'FtoS')
% %         inputs{1} = cellstr(structFile);
        inputs{1} = {fullfile(structPath,structf(1).name)};  % ref (fixed) -> struct
        file = functFilenames(1,:);
        %name=file.name;
        if RestingState_mode
            inputs{2} = {fullfile(alignFolder, ['mean' file.name(2:end)])}; % source (moved) -> meanfunc
        else
            inputs{2} = {fullfile(alignFolder, ['mean' file.name])};
        end
        display('\n CoregReslice FtoS \n');
        tmp_others= functf_fullpath';
%       tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'reg'; % prefix
    elseif strcmp(coregDirection,'StoF')
        inputs{1} = cellstr(refMeanFunc);  % ref (fixed) -> meanfunc
        file = functFilenames(1,:);
        inputs{2} = {fullfile(alignFolder, ['mean' file.name])}; % source (moved) -> func
        display('\n CoregResclice StoF \n');
        tmp_others=cellstr(spm_select('List',alignFolder,'rp*.txt'));
        tmp_others=cellfun(@(x) fullfile(alignFolder,x),tmp_others,'UniformOutput',false); % prepend dir
        inputs{3} = tmp_others; % others (all other funcs)
        inputs{4} = 'reg'; % prefix
    else
        error(['unknown coregDirection: ' coregDirection]);
    end
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.coregReslice = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% smooth

if  tasksTodo.smooth
    jobs = {fullfile(jobsParentPath,'jobs','smooth_job.mat')};
    inputs = cell(3, 1);
    % select all funcs
    tmp_data=cellstr(spm_select('List',alignFolder,'r*.nii'));
    tmp_data{end+1}=['mean' functFilenames(1,:)];
    tmp_data=cellfun(@(x) fullfile(alignFolder,x),tmp_data,'UniformOutput',false); % prepend dir
    inputs{1} = tmp_data;   % functional data
    inputs{2} = smoothFWHM; % smoothing kernel specs
    inputs{3} = ['s' num2str(smoothFWHM(1))];
    spm_jobman('serial', jobs, '', inputs{:});
    
    %---------------- update tasksTodo -------------------------
    copyfile(jobs{:}, jobFolder);
    tasksDone.smooth = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% Old segmentation and normalization (classic segmentation in SPM8)
% if tasksTodo.segment == OLDSEGMENT
%     structVol = spm_vol(structFile);
%     templVol = spm_vol(fullfile(spm('Dir'),'templates','T1.nii'));
%     
%     % -------------- segmentation --------------------------------
%     fprintf('%s','* Segmentation... ');
%     mkdir(segFolder);
%     defaults.segment.write.wrt_cor = 0;
%     
%     pm_segment(structVol,templVol,defaults.segment);
%     movefile(fullfile(structPath,'c*'), segFolder);
%     fprintf('%s\n','Segmentation done.');
%     
%     % -------------- normalization -----------------------------
%     fprintf('%s','* Normalisation... ');
%     mkdir(mniFolder);
%     
%     spm_normalise(templVol, structVol);
%     matfile = fullfile(structPath,[structFilename(1:end-4) '_sn.mat']);
%     spm_write_sn(structVol.fname, matfile);
%     
%     movefile(fullfile(structPath,'w*'), mniFolder);
%     movefile(fullfile(structPath,'*sn.mat'), mniFolder);
%     
%     fprintf('%s\n','Normalisation done.');
%     
%     % ------------ update tasksTodo -------------------------------
%     tasksDone.segment = OLDSEGMENT;
%     save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
%     
% end

%% Segmentation and normalization ("New segmentation" in SPM8)
if tasksTodo.segment == SEGMENT
    %structf = dir(fullfile(structPath,[prefix_struct,'*']));
    cd(alignFolder);
    structf = dir(fullfile(structPath,'*.nii'));
    % --------- do segmentation ---------------------------------
    fprintf('%s','* Segmentation... ');
    if ~exist(segFolder,'dir')
        mkdir(segFolder);
    end

    jobs = {fullfile(jobsParentPath,'jobs','segmentWithWarped_job.mat')};
    %jobs = {fullfile(jobsParentPath,'jobs','segmentWithWarped_job_linux.mat')};
    spm_jobman('serial', jobs, '', cellstr(fullfile(structPath, structf(1).name)));
    % arrange files
    movefile(fullfile(structPath,'i*'), segFolder);
    movefile(fullfile(structPath,'y*'), segFolder);
    movefile(fullfile(structPath,'c*'), segFolder);
    movefile(fullfile(structPath,'*seg8.mat'), segFolder);
    
    mkdir(mniFolder);
    movefile(fullfile(structPath,'w*'),mniFolder);
    if size(dir(fullfile(structPath,'mw*')),1)
        movefile(fullfile(structPath,'mw*'),mniFolder);
    end
    
    copyfile(jobs{:}, jobFolder);

    fprintf('%s\n','Segmentation done.');
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = SEGMENT;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% Segmentation + DARTEL
if tasksTodo.segment == SEGMENTDARTEL
    % --------- do segmentation ---------------------------------
    mkdir(segFolder);
    
    
    jobs = {fullfile(jobsParentPath,'jobs','segmentDartel_job.mat')};
    spm_jobman('serial', jobs, '', cellstr(structFile));
    
    % arrange files
    movefile(fullfile(structPath,'i*'), segFolder);
    movefile(fullfile(structPath,'y*'), segFolder);
    movefile(fullfile(structPath,'c*'), segFolder);
    movefile(fullfile(structPath,'rc*'), segFolder);
    movefile(fullfile(structPath,'*seg8.mat'), segFolder);
    
    copyfile(jobs{:}, jobFolder);
    
    % ------------ update tasksTodo -------------------------------
    tasksDone.segment = SEGMENTDARTEL;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end


%% IBASPM Labelling
if 	tasksTodo.label
    fprintf('%s','* Atlasing... ');
    cd(segFolder);
    list = dir('c*.nii');
    %list = spm_select('List',segFolder,'c*.nii');
    segmentFiles = cell(size(list,1),1);
    for i = 1:size(list,1)
        segmentFiles{i}= fullfile(segFolder, list(i).name);
    end
    
    if (tasksDone.segment == SEGMENT) || (tasksDone.segment == SEGMENTDARTEL), deform = fullfile(segFolder, ['iy_' structFilename]);
    elseif tasksDone.segment == OLDSEGMENT, deform = fullfile(normFolder, [structFilename(1:end-4) '_sn.mat']); end
    
    if atlasType(1)=='G' %Greicius
        %repeat the labeling for each atlas map
        for aaa=1:size(atlasFiles,2)
            new_Auto_Labelling(segmentFiles, atlasFiles{aaa}.fname, deform, atlasFolder, []);
        end
    end
    
    if atlasType(1)=='A' %AAL
        new_Auto_Labelling(segmentFiles, atlasFile, deform, atlasFolder, []);
    end
        
    
    fprintf('%s\n','Atlasing done.');
    
    % I comment it so I don't update the labeling job and I can redo it
    % with different atlases...Giulia
    %---------------- update tasksTodo -----------------------
    tasksDone.label = 1;
    save(fullfile(jobFolder, 'tasksDone.mat'), 'tasksDone');
end

%% end
Ttotal=etime(clock, Tstart);
disp(['** DONE PROCESSING. Total time: ' num2str(Ttotal/60,'%3.1f') ' min.']);




