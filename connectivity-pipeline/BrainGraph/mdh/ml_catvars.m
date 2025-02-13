function ml_catvars(var,dim)
% ML_CATVARS Concatinates variables in the calling workspace.
%   ML_CATVARS(var,dim)
%       Will concatinate the variables var0,var1,...,varn along the
%       dimension dim in the calling workspace. var0,var1,...,varn will be
%       cleared from the workspace and replaced with a single variable var.
%
%   See also ML_SPLITVAR.
%
%   Author:
%       Martin Larsson
%       March 2017

    i = 0;
    value = [];
    name = [var '0'];
    while evalin('caller',['exist(''' name ''',''var'')'])
        value = cat(dim,value,evalin('caller',name));
        evalin('caller',['clear(''' name ''')']);
        i = i+1;
        name = [var num2str(i)];
    end
    assignin('caller',var,value);
end

