function [f,logfile1,logfile2] = mdh_dsistudio_wb(f,dsi_root,varargin)
% MDH_DSISTUDIO Builds ODFs by running DSI STUDIO.

% HB [Feb 2018]
% - specifying dsi studio root
% - option to vary N_fibers
% - input/output file struture; more robust file namings.
% 
% Contributers: 
%    Martin Larsson
%    David Abramian
%    Hamid Behjat
% 2017-2018.

% Code revised by: Anjali Tarun:
% Changes made:
%     -- to use whole brain volume instead of WM only
%     -- obtains whole brain brain mask
% February 2018


cntr = {
    'N_fibers',5,... % maximum count of the resolving fibers
    };
argselectAssign(cntr);
argselectCheck(cntr,varargin);
argselectAssign(varargin);

f.source = fullfile(f.diffusion,'data.nii.gz');
f.out_src = fullfile(f.diffusion,'odfs.src.gz');
f.mask = fullfile(f.T1w,'brainmask_fs_125.nii');
f.odfs = fullfile(f.diffusion,'odfs.fib.gz');

d = fullfile(fileparts(f.out_src),['odfs.src.gz.odf8.f',...
    num2str(N_fibers),'rec.bal.012fy.rdi.gqi.1.25.fib.gz']);

if ~exist(f.out_src,'file')
    fprintf('..Running DSI Studio: obtaining .src file \n')
    command = [dsi_root,'/dsi_studio --action=src --source=',...
        f.source,' --output=',f.out_src];
    [sts,logfile1] = system(command);
    if ~sts
        error('error running dsi studio [constructing odfs]');
    end
    fprintf('..done.\n');
    
    fprintf('..Running DSI Studio: obtaining .fib file \n')
    command = [dsi_root,filesep,'dsi_studio --action=rec --source=',...
        f.out_src,' --thread=2 --mask=',f.mask,...
        ' --method=4 --param0=1.25 --odf_order=8 --num_fiber=',...
        num2str(N_fibers),' --check_btable=1 --record_odf=1'];
    [sts,logfile2] = system(command);
    if sts
%         out_string = 'output odfs';
%         I = strfind(logfile2,out_string) + length(out_string) + 1;
%         f.out_fib = logfile2(I:end-1);
%         copyfile(f.out_fib,f.odfs);
          outfile = dir(fullfile(f.diffusion,'odfs.src.gz.odf*'));
          copyfile(fullfile(f.diffusion,outfile(1).name),f.odfs);
    else
        error('error running dsi studio [constructing odfs]');
    end
    fprintf('..done.\n');
else
    fprintf('..ODF file already exists.\n')
end

