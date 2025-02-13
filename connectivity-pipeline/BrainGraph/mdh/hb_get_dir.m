function dirr = hb_get_dir(dirCatgory,varargin)
%
%
%
% Requires:
%  hb_get_tag.m
%  hb_get_name.m
%
% Hamid Behjat
% March-July 2017

cntr = {
    'baseDir',[],...           
    'subjectID',[],...  
    'subjectDir',[],... * for 'local_graphs'
    'opts',[],...
    'G_num',[],...
    'iG',[],...
    'N_regions',[],...
    'ribbonDir',[]};              
argselectAssign(cntr);
argselectCheck(cntr,varargin);
argselectAssign(varargin);

switch dirCatgory
    case 'hcp_t1w_mri'
        dirr = [baseDir,filesep,num2str(subjectID),filesep,'T1w'];%[12Feb2018]      ,filesep,num2str(subjectID),filesep,'mri'];
        
        % tempDir = [dataDirectory,filesep,num2str(tempIDs(iSubject)),filesep,'T1w',filesep,num2str(tempIDs(iSubject)),filesep,'mri'];
    case 'hcp_t1w_surf'
        dirr = [baseDir,filesep,num2str(subjectID),filesep,'T1w',filesep,num2str(subjectID),filesep,'surf'];
        
    case 'hcp_mni'
        dirr = [baseDir,filesep,num2str(subjectID),filesep,'MNINonLinear'];
    
    case 'bjt'
        dirr = [ribbonDir,filesep,'bjt'];
        
    case {'local_graphs','local_graphs_lr','local_graphs_nol'}
        % nol: non-overlapping
        
        dirr = [...
            subjectDir,...
            filesep,...
            hb_get_name('vol_clusters_long',opts,...
            'N_regions',N_regions,'justName',1)];
        
        if any(strcmp(dirCatgory,{'local_graphs','local_graphs_nol'}))
            dirr = [dirr,filesep,hb_get_tag('G_size',opts)];
        elseif strcmp(dirCatgory,'local_graphs_lr')
            dirr = [dirr,filesep,hb_get_tag('G_size_lr',opts)];
        end
        
        if any(strcmp(dirCatgory,{'local_graphs','local_graphs_lr'}))
            switch opts.localGraphConstructionScheme
                case 'nearestNeighbsAroundCenter'
                    dirr = [dirr,'nearestNeighbs_'];
                    
                case 'diffuseOutFromCenter'
                    dirr = [dirr,'diffuseOut_'];
                    
                otherwise
                    error('ooopps..')
            end
        elseif strcmp(dirCatgory,'local_graphs_nol')
            dirr = [dirr,'nonOverLapping_'];
        end
        
        dummy = hb_get_tag('A_graph_short',opts);
        dirr = [dirr,dummy(1:end-1)];
        
    case 'signals'
        [pathstr,dummy] = fileparts(hb_get_name('mat_local_graph_mini',opts,'subjectDir',subjectDir,'N_regions',N_regions,'iG',iG));
        dirr = [pathstr,filesep,dummy];
        
    case 'hcp_mni_fmri'
        dirr = [baseDir,filesep,num2str(subjectID),filesep,'MNINonLinear',filesep,'Results'];
        
    otherwise
        error('ooopps..')
        
end






























