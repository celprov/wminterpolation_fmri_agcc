% -------------------------------------------------------
%
%    pp_loadVolumes - modified from QNAP
%
%    Created:         Daniela Zoeller      (30.10.2015)
%    Based on: pp_loadVolumes.m (Jonas Richiardi and Giulia Preti)
%    Last modified:   Daniela Zoeller      (30.10.2015)
%
%    Medical Image Processing Lab
%    EPFL - UniGe
%
% ------------------------------------------------------
%
% Modification pp_loadVolumes.m:
%   30.10.2015: - add option newSegment with Dartel and no deformation 
%                   task and registration to MNI



% Helper script for preprocess.m - in charge of loading data
% v1.0 Manuel Wuethrich- initial release
% v1.0.1 JR - cleaner directory tests
% v1.0.2 JR - handle coregResliceFtoS
% v1.0.3 JR - handle coreg FtoS, smoothing
% v1.0.4 JR June 2011 - more flexible handling of coregistration
% v1.0.5 JR June 2011 - handle custom atlas file
% v1.0.6 JR Dec 2011 - checks specifically for files and dir types, not
%                       just their existence
% v1.0.7 JR Dec 2011 - corrects bug on tasksDone checking, more thorough
%                       type checking
%% initialize
RestingState_mode = 1;
pp12_defineConstants;
if RestingState_mode
    [volExt,chain,QCcoef, highpass_bool,coregDirection,refMeanFunc,...
        smoothFWHM,atlasFile, atlasType]=process_options(varargin,...
        'volExt','nii','procChain',{'QC','coregReslice'},...
        'QCcoef', 1.5,'highpass',0,'coregDirection','FtoS','refMeanFunc','',...
        'smoothFWHM',[6 6 6],'atlasFile',[],'atlasType','AAL');
else
    [volExt,chain,QCcoef, highpass_bool,coregDirection,refMeanFunc,...
        smoothFWHM,atlasFile, atlasType]=process_options(varargin,...
        'volExt','nii','procChain',{'reset','realign','QC','coregReslice'},...
        'QCcoef', 1.5,'highpass',0,'coregDirection','FtoS','refMeanFunc','',...
        'smoothFWHM',[6 6 6],'atlasFile',[],'atlasType','AAL');
end
%{'reset','realign','QC','coregister','smooth','segment','label'}
%% set folders
segFolder = fullfile(structPath,'Segmented');
dartelFolder = fullfile(structPath,'DARTEL');
mniFolder = fullfile(structPath,'MNI');
atlasFolder = fullfile(structPath,'Atlased');
alignFolder = fullfile(functPath,'realigned');
QCFolder = fullfile(functPath,'QC');

minLength = min(size(structPath,2),size(functPath,2));
structPath_exists=exist(structPath,'dir')==7;
functPath_exists=exist(functPath,'dir')==7;
if structPath_exists && ~functPath_exists, jobFolder = fullfile(structPath,'jobs');
elseif ~structPath_exists && functPath_exists, jobFolder = fullfile(functPath,'jobs');
elseif structPath_exists && functPath_exists
    if strcmp(structPath, functPath), jobFolder = fullfile(functPath,'jobs');
    else jobFolder = fullfile(fileparts(structPath(1:find(structPath(1:minLength) ~= functPath(1:minLength), 1 ))),'jobs');
    end
else error(['Neither functional' functPath ' nor structural path ' ...
        structPath  exist: dir(fullfile(alignFolder,['r*' volExt '$'])) ]);
end
if exist(jobFolder,'dir')~=7, mkdir(jobFolder); end;
% % CHECK ATLAS
% if isempty(atlasFile)
%     atlasFile = fullfile(spm('Dir'),'toolbox','IBASPM','atlas116.img');    
% end
% if exist(atlasFile,'file')~=2
%     error(['Atlas file ' atlasFile ' does not exist. Check your paths.']);
% end

%% save chain to the struct tasksTodo and check arguments
tasksTodo = struct('reset',0,'realign',0,'QC',0,'coregister',0,'coregReslice',0,...
    'smooth',0,'segment', 0,'label',0);
for procElem = chain
    switch procElem{1}
        case 'reset'
            tasksTodo.reset = 1;
        case 'realign'
            tasksTodo.realign = 1;
        case 'QC'
            tasksTodo.QC = 1;
        case 'coregister'
            tasksTodo.coregister = 1;
            if tasksTodo.coregReslice == 0
                tasksTodo.coregister = 1;
            else
                error('There must be only one coregistration step');
            end
        case 'coregReslice'
            tasksTodo.coregReslice = 1;
            if tasksTodo.coregister == 0
                tasksTodo.coregReslice = 1;
            else
                error('There must be only one coregistration step');
            end
        case 'smooth'
            tasksTodo.smooth=1;
        case 'oldSegment'
            if tasksTodo.segment == 0
                tasksTodo.segment = OLDSEGMENT;
            else
                error('There must not be more than one segmentation method.');
            end
        case 'segment'
            if tasksTodo.segment == 0
                tasksTodo.segment = SEGMENT;
            else
                error('There must not be more than one segmentation method.')
            end
        case 'segmentDartel'
            if tasksTodo.segment == 0
                tasksTodo.segment = SEGMENTDARTEL;
            else
                error('There must not be more than one segmentation method.')
            end
       case 'label'
            tasksTodo.label = 1;
        otherwise
            error(['processing command "' procElem{1} '" does not exist']);
    end
end

if exist(fullfile(jobFolder, 'tasksDone.mat'),'file')==2
    load(fullfile(jobFolder, 'tasksDone.mat'))
else
    tasksDone = struct('reset',0,'realign',0,'coregister',0,'coregReslice',0,...
        'smooth',0,'segment', 0,'label',0,'VERSION','2.1.2');
end

if tasksDone.realign, tasksTodo.realign = 0; end
%if tasksDone.reset, tasksTodo.reset = 0; end
if isfield(tasksDone, 'QC') && tasksDone.QC, tasksTodo.QC = 0; end
if tasksDone.coregister, tasksTodo.coregister = 0; end
if tasksDone.coregReslice, tasksDone.coregReslice = 0; end
if tasksDone.smooth, tasksTodo.smooth = 0; end
if (tasksTodo.segment && tasksDone.segment)
    if tasksDone.segment == tasksTodo.segment
        tasksTodo.segment = 0; %Daniela: comment to do same segmentation multiple times (jobs debugging)
    else
        warning('A different segmentation method has already been applied.')
    end
end
if tasksDone.label, tasksTodo.label = 0; end

%% arguments sanity check
if (length(volExt)~=3)
    error('Volume filename extension must be 3 characters long');
end

% if (tasksTodo.label && ~tasksTodo.segment && ~tasksDone.segment)
%     error('Before labelling volume must be segmented. Add "newSegment" or "classicSegment" to your processing chain.');
% end


% XXX 
%if (tasksTodo.coregister && ~tasksTodo.realign && ~tasksDone.realign)
%    error('Before coregistering functional volumes must be realigned. Add "realign" to your processing chain.');
%end
%if (tasksTodo.coregReslice && ~tasksTodo.realign && ~tasksDone.realign)
%    error('Before coregistering functional volumes must be realigned. Add "realign" to your processing chain.');
%end
% XXX ADD CHECK FOR SMOOTH - MUST BE REALIGNED AND COREGED FIRST
if (tasksTodo.QC && ~tasksTodo.realign && ~tasksDone.realign)
%    error('Before quality control functional volumes must be realigned. Add "realign" to your processing chain.');
end
if (tasksTodo.coregReslice && strcmp(coregDirection,'FtoF') && isempty(refMeanFunc))
    error('A reference functional image is required for FtoF coregReslice');
end

%% load atlas file/files
if tasksTodo.label
    
    % CHECK if ATLAS exist
    if exist(atlasFile,'file')~=2
        error(['Atlas file ' atlasFile ' does not exist. Check your paths.']);
    end
    
    % check if it is Greicius
    if atlasType(1)=='G'
        %enter inside folder and load all the Greicius maps in the cell
        %atlasFiles
        atlasFiles=[];
        lp=strfind(atlasFile,'/');
        lp=lp(end);
        atlasdir=atlasFile(1:lp-1);%directory with Greicius maps
        clear lp;
        [Greiciusmaps] = spm_select('List',atlasdir,['.*\.nii$']);
        if isempty(Greiciusmaps)
            [Greiciusmaps, foosubdirs] = spm_select('List',atlasdir,['.*\.img$']);;
        end
        for r=1:size(Greiciusmaps,1)
            atlasFiles{r}=spm_vol(fullfile(atlasdir,Greiciusmaps(r,:)));
        end
    end
end

%% check if files are compressed

%--------- structural files ---------------------------
niiDirs = dir(fullfile(structPath,['*.',volExt]));
gzDirs = dir(fullfile(structPath,'*.gz'));

if size(gzDirs,1) && ~size(niiDirs,1)
    gzStructFile = fullfile(structPath, gzDirs(1).name);
    gunzip(gzStructFile);
    
    vol = spm_vol(gzStructFile(1:end-3));
    vol.mat(:,4) = [130 -130 -60 1];
    %vol.fname=strcat(vol.fname(1:end-4),'f','.nii');
    spm_write_vol(vol, spm_read_vols(vol));
end

%--------- functional files ---------------------------
niiDirs = dir(fullfile(functPath,['*.',volExt]));
gzDirs = dir(fullfile(functPath,'*.gz'));

if size(gzDirs,1) && ~size(niiDirs,1)
    gzFunctFile = fullfile(functPath, gzDirs(1).name);
    gunzip(gzFunctFile);
end

%% check if 4d nii file
niiDirs = dir(fullfile(functPath,['*.',volExt]));

if size(niiDirs,1) == 1
    functFile = fullfile(functPath, niiDirs(1).name);
    vols = spm_vol(functFile);
    
    for i = 1:size(vols,1)
        vol = vols(i);
        [path name ext] = fileparts(vol.fname);
        numberStr = ['00' int2str(i)];
        numberStr = numberStr(size(numberStr,2)-2:end);
        
        %vol.mat(:,4) = [110 -110 0 1];
        
        vol.fname = fullfile(path, ['f' name numberStr ext]);
        vol.n = [1 1];
        fprintf(num2str(i))
        spm_write_vol(vol, spm_read_vols(vols(i)));
        
    end
end

%% select structural files
if tasksTodo.coregister || tasksTodo.coregReslice || tasksTodo.segment || tasksTodo.label
    if exist(structPath,'dir')~=7, error(['Structural path does not exist: ' path_s ]); end
    
%     structFilename = spm_select('List',structPath,'(^s|^t1|^ms)*.nii');
    cd(structPath);
    structFile = dir('rs*.nii');
    structFilename = structFile(1).name;
    %if size(structFilename,1)~=1, error('There should be exactly 1 structural volume, check code.'); end
    structFile = fullfile(structPath, structFilename);
end

% select functional files (Giulia:commented if next one is not commented)
if tasksTodo.realign || tasksTodo.coregister || tasksTodo.coregReslice ...
        || tasksTodo.QC || tasksTodo.smooth
    if exist(functPath,'dir')~=7, error(['Functional path does not exist: ' path_f ]); end
    
    if RestingState_mode
        cd(fullfile(functPath,'realigned')); 
        functFilenames = dir('r*.nii');
    else
        cd(functPath); 
        functFilenames= dir('*.nii');
    end
    
    nFiles=size(functFilenames,1);
    
    if nFiles==0 && (tasksTodo.realign || tasksTodo.coregister)
        error('No functional files selected');
    else
        disp(['Selected ' num2str(nFiles) ' functional volume files for processing']);
    end
    % prepend full dir;
    functFiles=cell(nFiles,1);
    for f=1:nFiles
        file = functFilenames(f,:);
        if RestingState_mode
            functFiles{f}=fullfile(functPath,'realigned',file.name);%resting
        else
            functFiles{f}=fullfile(functPath,file.name); %working
        end
    end
    

end

%% select ONLY mean functional file (Giulia realigned already done, I already have the mean AND DON'T WANNA LOAD FUNCTINAL FILES)
%     
%    functFilename = spm_select('List',strcat(functPath,'/realigned'),['^mean.*\.*']);
%     nFiles=size(functFilename,1);
%     if nFiles==0 && (tasksTodo.realign || tasksTodo.coregister)
%         error('No mean selected');
%     else
%         disp(['Selected mean functional image for processing']);
%     end
%     % prepend full dir;
%     meanFile=cell(1,1);
%     
%     meanFile=fullfile(strcat(functPath,'/realigned'),functFilename);
%   
%%
dispTasks12(tasksTodo, tasksDone);
