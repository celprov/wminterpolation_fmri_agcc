function [A,indices,G,A_a,A_ac,A_acp] = hb_adjacencymatrix_gm(mask,conn,connComp,varargin)
% HB_ADJACENCYMATRIX_GM Constructs gray matter graph A matrix.
%
% Inputs:
%   mask: 3D volume, or its full address (/Users/.../xxx.nii).
%   conn: 3D neighbourhood: 6,18 or 26. 
%   connComp: Number of connected components in structure.  
%   
% Outputs:
%   A: adjacency matrix
%   indices: corresponding indices in 'mask'
%   Aa:
%   Aac:
%   Aacp:
% 
% Contributers: 
%    Hamid Behjat
%    Martin Larsson
% March 2017 

cntr = {
    'weight','no',...
    'pialFiles',[],...
    'parallelize',0,...
    'sourceFile',[],...
    'G',[]}; 
argselectAssign(cntr);
argselectCheck(cntr,varargin);
argselectAssign(varargin);

fprintf('\n.Constructing GM graph adjacency matrix..\n')
if ~isempty(pialFiles) 
    if ischar(mask)
            sourceFile_h = spm_vol(mask);
    elseif ~isempty(sourceFile)
        sourceFile_h = spm_vol(sourceFile);
    else
        error('source file needs to be specified to extract header.')
    end
end

% compute adjacency mtrix (A)
fprintf('..Determining adjacencies..\n')
[A,indices,mask] = hb_compute_adjacency(mask,conn,...
    'weight',weight,'sourceFile',sourceFile,'sS',0);
fprintf('..done.\n')

A_a = A;  % A adjacency
indices_a = indices; 

% Clean A.
[A,d,mask] = ml_clean_adjacency(A,connComp,mask);
%fprintf([num2str(numel(indices)-numel(d)),' nodes removed.\n'])

A_ac = A; % A adjacency-cleaned

indices = indices(d);

% sanity check
if setdiff(indices,find(mask))
    error('something is fishy..')
end

if ~isempty(pialFiles)
    pials = ml_batch_gifti(pialFiles,sourceFile_h.mat);
    A_tmp = A;
    
    % Prune A.
    fprintf('..Pruning adjacency matrix..\n')
    if conn==6
        A = ml_prune_adjacency(A,mask,pials,...
            'parallelize',parallelize,'conn6',1);
    else
        A = ml_prune_adjacency(A,mask,pials,...
            'parallelize',parallelize);
    end
    fprintf([' ',num2str((numel(find(A_tmp))-numel(find(A)))/2),...
        ' edges removed.\n'])
    fprintf('..done.\n')

    A_acp = A; % A adjacency-cleaned-pruned
    
    % Clean A.
    fprintf('..Cleaning adjacency matrix..\n')
    [A,d,mask] = ml_clean_adjacency(A,connComp,mask);
    fprintf([' ',num2str(numel(indices)-numel(d)),' nodes removed.\n'])
    fprintf('..done.\n')

    indices = indices(d);
    
    % sanity check
    if setdiff(indices,find(mask))
        error('something is fishy..')
    end
end

if ~isempty(G)
    G.N_gm = size(A,1);
    G.A_gm = A;
    G.indices_gm = indices;
    G.A_gm_nonpruned = A_a;
    G.indices_gm_nonpruned = indices_a;
end

A_gm = A; %#ok<NASGU>
indices_gm = indices; %#ok<NASGU>
save(fullfile(G.f.graph,'A_gm.mat'),'A_gm');
save(fullfile(G.f.graph,'indices_gm.mat'),'indices_gm');
fprintf('done.\n')
