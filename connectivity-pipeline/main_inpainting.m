%% Methodological Pipeline for the Manuscript

%     Structural mediation of human brain activity revealed by white-matter 
%     interpolation of fMRI
% 
%     Created by: Anjali Tarun
%     Date created: 21 February 2018
%     Modified: 7 August 2019
%              -- modified and organized for manuscript submission
%     DESCRIPTION:
% 
%         This is the main manuscript for performing white matter
%         inpainting from fMRI data using a voxel-wise brain connectome.
%         The pipeline works as follows:
% 
%                 1. The voxel-wise brain grid is first constructed. This
%                 is done by first extracting the ODF from the diffusion
%                 data. After which, the brain graph is built and the 
%                 Laplacian is computed. 

%                 2. Graph signal recovery is performed by minimizing a
%                 cost function that is described in the manuscript. 
%                 The cost function balances between (1) minimizing the RSS
%                 and retaining the original GM signals, and (2) imposing
%                 smoothness with respect to the structure of the graph. 
% 
% 
%     The pipeline is adjusted to take HCP folder as input and outputs 
%     the brain graph adjacency matrix, Laplacian, and the corresponding
%     interpolated volumes for each subject and session.
% 
%     Input:  HCP folder containing subject IDs
%             The HCP folder should contain (1)preprocessed diffusion data (2)
%             processed fMRI volumes (3) Extended structural volumes
%     Output: Voxel-brain graph and interpolated volumes
% 
%     The struct 'param' contains all the parameter information and paths
    

%% Specifies code paths and obtain data paths
clear all; close all;
%warning('off');

dataset = 'RestingState';

% celinepath = 'J:\Anjali_Diffusion_Pipeline\Celine';
celinepath = '/media/miplab-nas2/Data/Anjali_Diffusion_Pipeline/Celine';
param.spm = fullfile(celinepath,'spm12');
param.DSIstudio = fullfile(celinepath,'connectivity-pipeline','dsi_studio_64');
param.HCPDatapath = '/media/miplab-nas2/Data/AgCC/Anjali/BrainGraph_results';
param.HCPDatapath2 = fullfile(celinepath,'data',dataset);

addpath(genpath(param.spm))
addpath(genpath(pwd));

%% Find the subjects folder
dirs = dir(fullfile(param.HCPDatapath2,'s*'));
dirs = struct2cell(dirs);
dirs = dirs(1,:)';
param.Subjects = dirs;

Mask = zeros(1,length(dirs));
for i = 1:length(dirs)
    subject = param.Subjects{i};
    if exist(fullfile(param.HCPDatapath,subject, 'T1w','new_brainmask.nii'),'file')
        Mask(i) = 1;
    end
end

dirs = dirs(logical(Mask));
param.Subjects = dirs;        

for iS = 1:length(param.Subjects)
        tic
        param.subject = param.Subjects{iS};
        disp(['Constructing brain graph for subject..',param.subject])
       
        %% Construct the voxel-level brain grid from diffusion MRI

        % Calls for the parameters required for constructing the brain grid
        Inputs_BrainGrid

        % Runs the construction of the brain graph
        param = RunBrainGraph(param);

        %% Initialize the graph signal inpainting

        % First checks if the brain graph is already constructed
        if ~isfield(param.ODF.WB,'A')
            disp('Brain grid not available. Run construction of brain graph first!')
        end
        
        % Load parameters for functional volume masks
        Inputs_fMRI_Mask;
        param.SaveInpainting = fullfile(param.SaveInpainting, param.subject);
        if ~exist(param.SaveInpainting,'dir')
            mkdir(param.SaveInpainting)
        end
        
       % Initialize L2 norm problem    
        A_n = param.ODF.WB.A_n;
        param.indices_wb = indices_wb;
        n = size(A_n,1);

        % Construct the indicator matrix
        B = sparse(BC,BC,ones(1,length(BC)),n,n);

        BB = B'*B;

        L =  speye(size(A_n,1)) - A_n;

        % Computes for the parameter lambda to use
        lambda = CheckParameterLambda(param,BB,L,BC);
        
        %% Actual run of the signal inpainting after initialization

            % Load inputs for functional MRI
            % Pipeline assumes that the fMRI volumes have already been
            % preprocessed -- with prefix s6 (smooothed, Gaussian FWHM 6 mm)
          
            Inputs_fMRI
            
            disp(['Runing signal recovery for subject..',param.subject])
            
            %if ~exist(fullfile(param.SaveInpainting, ['GSR_Lambda',num2str(lambda),'_',param.subject,'.mat']),'file')
                
                disp(['Running for lambda..',num2str(lambda)])
                
                L_mat = BB+lambda*L;
                
                L1 = ichol(L_mat);
%                 alpha = max(sum(abs(L_mat),2)./diag(L_mat))-2;
%                 L1 = ichol(L_mat, struct('type','ict','droptol',1e-3,'diagcomp',alpha));

                V_new = zeros(NumScans,size(L1,1));

                % Main loop for signal inpainting
                % This may take a while -- alternatively one can run this
                % in parallel using parfor
                
                parfor i = 1:NumScans
                    disp(['Running for iteration..',num2str(i)])
                    xsignal = zeros(size(L1,1),1);
                    
                    xsignal(BC) = V(i,BC);

                    [V_new(i,:),~,~,~,~] = pcg(L_mat,xsignal,1e-8,100,L1,L1');
                    
                end
                
                V = V_new;
                clear V_new
                save(fullfile(param.SaveInpainting, ['GSR_Lambda',num2str(lambda),'_',param.subject,'.mat']),'V','-v7.3')
%                 clear V;     
        %end
        
        %% Visualize inpainted signal
        % save volumes here
%         path = fullfile(param.SaveInpainting, param.subject);
%         if ~exist(path, 'dir')
%             mkdir(path)
%         end
        
        for iter = 1:NumScans
            
            vol = zeros(fHeader.dim);
            vol(indices_wb) = V(iter,:);
            nfile = fullfile(param.SaveInpainting,['Inpainting',sprintf('_%.3d',iter),'.nii']);
            cbiWriteNifti(nfile,vol,hdr,'float32');
        end
        createfHeader(dataset);
      toc
end