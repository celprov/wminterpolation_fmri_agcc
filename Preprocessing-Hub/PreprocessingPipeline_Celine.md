# Copying the original files to your own folder
1. Run `arrange_data.m` and adapt the paths to your own needs

OR

1. Run `OrganizeNifti.m` to reorganize the emplacement of the nifti files


# Preprocessing

1. Run `pp12Main.m` using the options {'reset','realign','QC','coregReslice','segment'}  
	with parameters {'QCcoef',1.5,'coregDirection','FtoS','refMeanFunc',...
	'atlasFile',AALfile,'atlasType','AAL','highpass',0}

	Attention to change the RestingState_mode in pp12_loadVolumes 
	and preprocess12 as well !

	If you have errors, they most probably come from the fact that the files are not 
	found, either because they don't exist, either because the paths are incorrect.
	This applies even though the errors look more complicated than a file is not found.

	For certain subject coregistration called by pp12Main.m may fail, i.e coregistration to the T1 image
	that itself has been coregistered onto the diffusion data. In this case, run `different_coregistration.m` 
	to redo the coregistration using coregistration onto the diffusion data directly.

2. Run `PreprocROITimeCourse_AT.m` in the folder `NuisanceRegression` to perform nuisance regression

	Attention this needs to be run with SPM8 !