function [A,indices,mask] = hb_compute_adjacency(mask,conn,varargin)
% Computes the adjacency matrix for a given 3D (i.e. volume) mask.
%
% mask: input mask to compute its adjacency matrix (A).
%           A will be a square, symmetric matrix. Its dimension is equal
%           to the noumber of non-zero elements in the input mask.
%
%           The input mask should be a volume. For an image, append
%           the image by two zero plane images to convert it to a volume.
%
%           * 'mask' can also be the full adress to a .nii file, which will 
%           then be loaded. 
%
% conn: neighbourhood connectivity to be considered when compuing
%          the edges
%
% - Optional inputs -
% 'weight': weighted or non-weighted A. Default: 'no'.
%              For an example weigted option us: 'yes_power_additive'
%
% 'pow'   : the power used in defining weigting for 'yes_power_additive'
%             weigting option.
%
% Output:
% A        : Adjacency matrix in sparse format; i.e. only non zeros elements are saved.
%             To view the full matrix, use the 'full' command.
%
% ind_rmd  : Removed indices after pruning. if no pruning done, [] returned. 
%
%
% Examples:
%
% Ex.1 Computes non-weighted A, with neighbourhood connectivity 26:
%
%       A=compute_adjacency(mask,26);
%
% Ex.2 Computes non-weighted A, with neighbourhood connectivity 6:
%
%       A=compute_adjacency(mask,6);
%
% Ex.3 The same as in #1:
%
%       A=compute_adjacency(mask,26,'weight','no');
%
% Ex.4 The same as in #1, but with weighted edges, using a power law:
%
%       A=compute_adjacency(mask,26,'weight','yes_power_additive','pow',3);
%
% Ex.5 Define a mask and compute its A matrix:
%
%       mask = zeros(3,3,3);
%       mask(:,:,1) = [1 1 1; 0 1 0; 0 0    0];
%       mask(:,:,2) = [0 0 0; 0 1 0; 0 0    0];
%       mask(:,:,3) = [0 0 0; 0 1 0; 1 0.5 0.5];
%       A_v1 = compute_adjacency(mask,26);
%       A_v2 = compute_adjacency(mask,6);
%       A_v3 = compute_adjacency(mask,26,'weight','yes_power_additive');
%       full(A_v1)
%       full(A_v2)
%       full(A_v3)
%
% Ex.6 Convert a 2D mask to 3D and compute its A matrix:
%
%       Im = [1 1 1 1 1; 0 0 0 1 0; 0 0 1 0 0; 0 1 0 0 0; 1 1 0.5 0.5 0.5]
%       mask = zeros(5,5,3);
%       mask(:,:,1) = zeros(5,5);
%       mask(:,:,2) = Im;
%       mask(:,:,3) = zeros(5,5);
%       A_v1 = compute_adjacency(mask,26);
%       A_v2 = compute_adjacency(mask,6);
%       A_v3 = compute_adjacency(mask,26,'weight','yes_power_additive','pow',5);
%       full(A_v1)
%       full(A_v2)
%       full(A_v3)
%--------------------------------------------
% 	Author: Hamid Behjat
%
%   Biomedical Signal Processing Group,
%   Dept. of Biomedical Engineering,
%   Lund University, Sweden
%
%   Oct 2016 / Feb 2018
%
% Acknowledgment:
% The initial version of part of this function was developed by
% Elena Najdenovska and Nora Leonardi, 2011.
% ************************************************

control_params = {
    'weight','no',...
    'pow',5,...
    's',[],...
    'distMetric','Eucl',...
    'sw',1,...
    'indGM',[],...
    'sourceFile',[],... % file associated to the input mask
    'sS',1};         
argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if ischar(mask)
    if sS, fprintf(['- source file: ',mask,'\n']); end
    mask = spm_read_vols(spm_vol(mask));
elseif ~isempty(sourceFile)
    if sS, fprintf(['- source file: ',sourceFile,'\n']); end
end

dim = size(mask);

N = numel(mask); 

indices = find(mask);

[alli,allj,allk] = ind2sub(dim,indices);

switch conn
    case 6
        nN = 3;
    case 18
        nN = 9;
    case 26
        nN = 13;
    otherwise
        error('undefined 3D connectivity neighbourhood!')
end

ci = repmat(alli,nN,1);
cj = repmat(allj,nN,1);
ck = repmat(allk,nN,1);

switch conn
    case 6
        ni = [alli  ; alli+1; alli  ];
        nj = [allj+1; allj  ; allj  ];
        nk = [allk  ; allk  ; allk+1];
    case 18
        ni = [alli  ;alli+1;alli+1;alli+1;alli  ;alli  ;alli  ;alli+1;alli+1];
        nj = [allj+1;allj  ;allj-1;allj+1;allj  ;allj+1;allj+1;allj  ;allj  ];
        nk = [allk  ;allk  ;allk  ;allk  ;allk+1;allk-1;allk+1;allk-1;allk+1];
    case 26
        ni = [alli;alli+1;alli+1;alli+1;alli;alli;alli;alli+1;alli+1;alli+1;alli+1;alli+1;alli+1];
        nj = [allj+1;allj;allj-1;allj+1;allj;allj+1;allj+1;allj;allj;allj-1;allj+1;allj+1;allj-1];
        nk = [allk;allk;allk;allk;allk+1;allk-1;allk+1;allk-1;allk+1;allk-1;allk-1;allk+1;allk+1];
end

maskZ = cat(2,mask,zeros(dim(1),1,dim(3)));
maskZ = cat(1,maskZ,zeros(1,dim(2)+1,dim(3)));
maskZ = cat(3,maskZ,zeros(dim(1)+1,dim(2)+1,1));

valid = (ni>=1 & ni<=dim(1) & nj>=1 & nj<=dim(2) & nk>=1 & nk<=dim(3));
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

tt = sub2ind(size(maskZ),ni,nj,nk);
ee = maskZ(tt);
valid = logical(ee);
ni = ni(valid);
nj = nj(valid);
nk = nk(valid);
ci = ci(valid);
cj = cj(valid);
ck = ck(valid);

cInd = sub2ind(dim,ci,cj,ck);
nInd = sub2ind(dim,ni,nj,nk);

switch weight
    
    case 'no'
        % non-weighted adjacency matrix, as used in Behjat et al., 2015 & 2016.
        if sS, fprintf('- weigting: none.\n'); end
        
        A = sparse([cInd,nInd],[nInd,cInd],ones(1,2*numel(ni)),N,N);
        
    case 'yes_power_additive'
        % Weighted adjacency matrix, as used in Behjat et al., 2013 & 2014.
        % The GM probability of nodes connected through an edge are used
        % to defined a weight for the edge.
        
        pro = 1;
        p = (mask(cInd).*mask(nInd)*pro).^pow;
        A = sparse([cInd,nInd],[nInd,cInd],[p,p],N,N);
        
    case 'yes_power_subtractive'
        
        pro = 1;
        p = (mask(cInd).*mask(nInd)*pro).^(-pow);
        A = sparse([cInd,nInd],[nInd,cInd],[p,p],N,N);
        
    case 'yes_Nora'
        
        d = mask(cInd)-mask(nInd);
        if isempty(s)
            s = mean(abs(d))*sw;
        end
        
        if strcmp(distMetric,'Gaus')
            sim = exp(-d.^2/(2*s^2));
        else
            sim = 1./(1+d.^2/(2*s^2));
        end
        A = sparse([cInd,nInd],[nInd,cInd],[sim,sim],N,N);
        
        
    case 'yes_distance'
        if sS, fprintf('- weigting: inverse Euclidean ditance.\n'); end

        [aha1, aha2, aha3] = ind2sub(dim, cInd);
        [nababa1, nababa2, nababa3] = ind2sub(dim, nInd);
        d = sqrt(sum(([aha1-nababa1,aha2-nababa2,aha3-nababa3]).^2,2));
        
        if isempty(s)
            s = mean(abs(d))*sw;
        end
        
        if strcmp(distMetric,'Gaus')
            sim = exp(-d.^2/(2*s^2));    % gaussian similarity
        elseif strcmp(distMetric,'adjEucl')
            sim = 1./(1+d.^2/(2*s^2));   % adjusted Euclidean similarity
        elseif strcmp(distMetric,'Eucl')
            sim = 1./d;                  % Euclidean similarity
        end
        A = sparse([cInd,nInd],[nInd,cInd],[sim,sim],N,N);
        
    case 'yesIndexed'  % Oct 2016
        
        d1 = mask(cInd)-mask(nInd);
        d2 = find(d1);
        p = zeros(size(d1));
        p(d2) = -1; %#ok<FNDSB>
        A=sparse([cInd,nInd],[nInd,cInd],[p,p],N,N);
end
A = remove_empty_rows_cols(A,indices,'maintain_inds');
end

function [X,inds,inds_rmd,c] = remove_empty_rows_cols(X,inds,removeType)

switch removeType
    case 'maintain_inds'
        c=find(~sum(X,1));
        c=c(~ismember(c,inds));
        X(:,c)=[];
        X(c,:)=[];
    case 'update_inds'
        c=find(~sum(X,1));
        X(:,c)=[];
        X(c,:)=[];
        inds_rmd = inds(c);
        inds(c)=[];
end
end

