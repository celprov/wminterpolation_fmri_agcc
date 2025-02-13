function [outputFileName,sts,Nr] = hb_make_single_connected(inputFileName,conn,varargin)
%  Checks whether there is only one connected structure in the volume. 
%  If the structure in the input volume is not single connected, a new 
%  volume will be save in the same location as the original volume which 
%  is connected. 
%  
%  
% Inputs:
%   inputFileName: Name of the source .nii file to be checked.
%                  The name should be full, i.e. should include full
%                  diectory address, file name and .nii at the end.
%
%   conn:          Desired connectivity. conn may have the following 
%                  scalar values:  
%
%                  4     two-dimensional four-connected neighborhood
%                  8     two-dimensional eight-connected neighborhood
%                  6     three-dimensional six-connected neighborhood
%                  18    three-dimensional 18-connected neighborhood
%                  26    three-dimensional 26-connected neighborhood 
%
%   varargin{1}:   Name of the .nii file to be saved.
%                  The name should be full, i.e. should include full
%                  diectory address, file name and .nii at the end.
%
% Outputs: 
%   outputFileName: Name of the saved .nii file. The name is full, i.e. 
%                   includes full diectory address, file name and .nii 
%                   at the end.
%
%   sts:            1 if there were voxels removed and the volume saved 
%                   is diferent from the original volume, otherwise 0.  
%
%   Nr:             Number of voxels removed.  
%
% Hamid Behjat
% Jan. 2017 / Feb 2018

[p,n,e] = fileparts(inputFileName);
control_params = {
    'outputFileName',fullfile(p,[num2str(conn),'connected_',n,e]),...
    'alwaysWriteConnected',0,'sS',1};

argselectAssign(control_params);
argselectCheck(control_params,varargin);
argselectAssign(varargin);

if sS, fprintf('\n* Making the structure a single-connected one..\n'); end

vol = spm_read_vols(spm_vol(inputFileName));
N_v0 = numel(find(vol(:)));

% remove possible small isolated objects
objSize = floor(numel(find(vol(:)))/3); % a large enough object size
vol = bwareaopen(vol,objSize,conn);
N_v1 = numel(find(vol(:)));
Nr = N_v0-N_v1;

if sS
    fprintf([' # of voxels removed to make struture single-connected: ',...
        num2str(Nr),'\n']);
end

if Nr~=0 || alwaysWriteConnected==1
    sts=1;
    V = spm_vol(inputFileName);
    V.fname = outputFileName;
    V.private.dat.fname = V.fname;
    V.descrip = [V.descrip,' -- made connected (conn: ',num2str(conn),')'];
    if sS, fprintf('- writing volume on disc: \n'); end
    spm_write_vol(V,vol);
    if sS, disp(V), end
else
    sts = 0;
end
if sS, fprintf('done.\n'); end
