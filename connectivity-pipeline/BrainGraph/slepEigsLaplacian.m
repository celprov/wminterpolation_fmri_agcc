function [Utr,S1]=slepEigsLaplacian(A,CONST_W,opts)

D = speye(size(A,1));

[Utr,S1]=eigs(@Binline,size(A,1),max(CONST_W),'sa',opts);
% [Utr,S1]=eigs(@Binline,size(A,1),max(CONST_W),'sr',opts);



function y=Binline(x)
    y=D*x-A*x;
end

end
