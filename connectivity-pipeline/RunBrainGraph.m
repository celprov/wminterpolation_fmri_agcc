function param = RunBrainGraph(param)


    % Date and time when the routines are called
    param.date = strrep(strrep(datestr(datetime('now')),' ','_'),':','_');
    
    param.constW = param.c_bandwidth;
    % Original parameters
    param_init = param;
        
    disp(['Extracting whole brain graph using ODF-based design for ..',param.subject])


    % Calls the mdh routine to construct the graph and return adjacency matrix
    
    if ~exist(fullfile(param.HCPDatapath,param.subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','A_wb.mat'))
        [param.ODF.WB,param.ODF.G] = main_braingraph(param);
        A_wb = param.ODF.WB.A;
        indices_wb = param.ODF.WB.indices_wb;
        
        save(fullfile(param.HCPDatapath,param.subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','A_wb.mat'),'A_wb','-v7.3')
        save(fullfile(param.HCPDatapath,param.subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','indices_wb.mat'),'indices_wb','-v7.3')
        
    else
        disp('Adjacency matrix already extracted..')
        load(fullfile(param.HCPDatapath,param.subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','A_wb.mat'))
        load(fullfile(param.HCPDatapath,param.subject,'ODF_Neigh_3_ODFPower_40','Spectrum_WB_improved','indices_wb.mat'))
        param.ODF.WB.A = A_wb;
        param.ODF.G.indices_wb = indices_wb;
    end
    
    [param.ODF.WB.A_n,param.ODF.WB.D]=slepNormalize(param.ODF.WB.A,param.normalize, param.normalize_type);

    % Computing the eigendecomposition of the Laplacian

    if param.percent
        param.constW = round(param.bandwidth*param.G.N_wb);
    else
        param.constW = param.c_bandwidth;
    end

    if param.Decomposition
        disp('Taking the eigenvalues...')

        [param.ODF.WB.Utr,param.ODF.WB.S]=slepEigsLaplacian(param.ODF.WB.A_n,param.constW,param.opts);

        % Saves the eigenvalues and eigenvectors

        SaveValues(param)

        SaveToNifti(param, param.ODF.WB.Utr,'Spectrum_WB')

    end
end