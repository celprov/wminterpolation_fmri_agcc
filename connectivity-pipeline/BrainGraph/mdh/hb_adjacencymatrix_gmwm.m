function [A,indices,G] = hb_adjacencymatrix_gmwm(G)
%HB_ADJACENCYMATRIX_GMWM Constructs GM+WM graph A matrix.
%   
%
%   Author:
%       Hamid Behjat
%       Feb 2018.

fprintf('\n.Constructing GM+WM graph adjacency matrix..\n')

[A,indices] = hb_compute_adjacency(G.f.gmwm_125,26,'sS',0);

[~,I_gm] = intersect(indices,G.indices_gm);
[~,I_wm] = intersect(indices,G.indices_wm);

d1 = setdiff(indices,G.indices_gm);
d2 = setdiff(d1,G.indices_wm);
[~,I_extra] = intersect(indices,d2);

% sanity checks
if length(I_gm)~= G.N_gm, error('fishy..'); end
if length(I_wm)~= G.N_wm, error('fishy..'); end
if (G.N_gm+G.N_wm+length(I_extra))~=length(indices), error('fishy..'); end

G.N_gmwm = G.N_gm+G.N_wm;

d1 = nnz(A); %#ok<NASGU>
A(I_extra,:) = 0;
A(:,I_extra) = 0;

% Find GM-WM boundary vertices.
d = A;
d(I_gm,I_gm) = 0;
d(I_wm,I_wm) = 0;
[d,~] = ind2sub(size(d),find(d));
d = unique(d);
d = indices(d);
[~,I_gm_boundary] = intersect(d,G.indices_gm);
[~,I_wm_boundary] = intersect(d,G.indices_wm);
if ~isempty(intersect(I_gm_boundary,I_wm_boundary)), error('fishy...'); end
G.indices_gm_boundary = d(I_gm_boundary);
G.indices_wm_boundary = d(I_wm_boundary);

% Restore GM edges.
d2 = nnz(A); %#ok<NASGU>
A(I_gm,I_gm) = G.A_gm;
d3 = nnz(A); %#ok<NASGU>

% Restore WM edges.
A(I_wm,I_wm) = G.A_wm;
d4 = nnz(A); %#ok<NASGU>

% Clean-up A.
c = find(~sum(A,1));
A(:,c) = [];
A(c,:) = [];
indices(c)=[];

% Store.
G.A_gmwm = A;
G.indices_gmwm = indices;
fprintf('.done.\n')

% Stuff.
d5 = nnz(A);
d6 = nnz(G.A_gm)/d5;
d7 = nnz(G.A_wm)/d5;
d8 = 1-d6-d7;

% Summary.
disp(' ');
disp(' ');
disp('*** GM+WM graph facts ***')
disp(['percentage of nodes being GM: ',num2str(100*G.N_gm/G.N_gmwm),'%']);
disp(['percentage of nodes being WM: ',num2str(100*G.N_wm/G.N_gmwm),'%']);
disp(' ');
disp(['percentage of GM nodes at GM-WM boundary: ',...
    num2str(100*length(G.indices_gm_boundary)/G.N_gm),'%']);
disp(['percentage of WM nodes at GM-WM boundary: ',...
    num2str(100*length(G.indices_wm_boundary)/G.N_wm),'%']);
disp(' ');
disp(['percentage of edges being intra GM   : ',num2str(100*d6),'%']);
disp(['percentage of edges being intra WM   : ',num2str(100*d7),'%']);
disp(['percentage of edges being inter GM-WM: ',num2str(100*d8),'%']);
disp(' ');
disp(['percentage of edges pruned: ',num2str(round(1e4*((d2-d3)/d2))*1e-4*100),'%']);

