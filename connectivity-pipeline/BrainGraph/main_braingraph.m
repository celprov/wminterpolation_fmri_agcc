%%    Main code for constructing the voxel-level brain grid

%     Created by: Anjali Tarun
%     Date created: 21 February 2018
%     Based from the original code of Hamid Behjat, Martin Larsson and
%     David Abramian


%     Changes made in the original ODF design of Itturia Medina and HB's code:
%     -- to use whole brain volume instead of WM only
%     -- incorporated structural boost using QA anisotropy
%     -- removes outlier brain voxels to aid in convergence of the
%           eigendecomposition of the Laplacian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [WB, G] = main_braingraph(param)

    hcp_root = param.HCPDatapath; % folder in which your HCP subjects are stored in

    dsi_root = param.DSIstudio; % root to dsi_studio's executable

    subjID   = param.subject;
      
    % Settings. 
    odfPow       = param.ODF.odfPow; % ODF power.
    neighborhood = param.neighborhood; % 3 or 5 [This is for white matter; gray matter always 3]
    N_fibers     = param.N_fibers; % dsi studio's and M&L's default: 5

    % Preprocess.
    G = mdh_preprocess(subjID,hcp_root);

    % Whole brain graph.
    [A_wb,indices_wb,G,O] = mdh_adjacencymatrix_wb(G,odfPow,neighborhood,...
        dsi_root,param,'N_fibers',N_fibers);

%     save(fullfile(G.f.graph,'G_wb.mat'),'G');

    WB.indices = indices_wb;
    WB.A = A_wb;
    WB.O = O;
    
   
end
