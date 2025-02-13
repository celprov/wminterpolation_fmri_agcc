
% Contains the input parameters for the construction of the brain grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if you want to get the eigenmodes of the brain graph put 1
param.Decomposition = 0;

param.ODF.odfPow = 40; % ODF power.
    
%     (!!!) It is very important to verify the neighborhood you wan to run
%     as it will dictate all the processes

param.ODF.neighborhood = 3; % 3 or 5 
param.ODF.N_fibers = 3; % dsi studio's and M&L's default: 5
% Tune FA/QA
param.alpha = 2;

% Choose the parameters for the eigendecomposition of the Laplacian to the 
% obtain eigenmodes..
% whether to have a percent bandwidth or constant
param.percent = 0; % 1 or 0
param.bandwidth = 1e-3; % if percent, we choose the percentage (*times the number of nodes)
param.c_bandwidth = 1000; % if constant, we choose the number of eigenmodes to extract
param.opts.issym=1;
param.opts.isreal=1;
param.opts.maxit=2500;
param.opts.disp=1;
param.normalize = 1;
param.normalize_type = 1; % The Laplacian is symmetrically normalized



% Specifies where to save all results of the brain grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

param.ODF.title = ['ODF_Neigh_',num2str(param.ODF.neighborhood),'_ODFPower_',num2str(param.ODF.odfPow)];
% Specifies where to save all results
param.ParentDirectory = fullfile(param.HCPDatapath2, 'BrainGraph_results');
param.ODF.SaveDirectory = fullfile(param.ParentDirectory,param.subject,param.ODF.title);
param.structural = fullfile(param.HCPDatapath2,param.subject, 'struct');
param.structseg = fullfile(param.structural,'Segmented');
if strcmp(param.subject, 's011')
    dirs = dir(fullfile(param.structural,'t*'));
else
   dirs = dir(fullfile(param.structural,'rs*')); 
end
param.T1file = dirs(1).name;
