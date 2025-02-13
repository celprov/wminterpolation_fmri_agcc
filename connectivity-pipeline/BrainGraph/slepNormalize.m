% Construct normalized (sparse) adjacency matrix
% 
% returns
% A - normalized adjacency matrix
% D - degree vector
%
% Dimitri Van De Ville, Sep 2014

% Modified by Anjali Tarun, include random walk normalization, April 3 2016

function [A,D]=slepNormalize(A,CONST_NORMALIZE,CONST_NORMALIZE_type)

if ~issparse(A)
    warning('Adjacency matrix A is not sparse.');
end

msize=size(A,1);

% normalize adjacency matrix
if CONST_NORMALIZE
    D=sum(A);
    
    switch CONST_NORMALIZE_type
        case 1
            Dn=1./sqrt(D);
            Dn=spdiags(Dn.',0,spalloc(msize,msize,msize));
            A=Dn*A*Dn;
        case 2
            Dn=1./D;
            Dn=spdiags(Dn.',0,spalloc(msize,msize,msize));
            A=Dn*A;
    end
end

% construct D
D=sum(A);
D=spdiags(D.',0,spalloc(msize,msize,msize));
