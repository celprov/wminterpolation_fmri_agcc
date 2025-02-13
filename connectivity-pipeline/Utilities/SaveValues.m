function [] = SaveValues(param)
     
    
    if ~exist(fullfile(param.SaveDirectory), 'dir')
         mkdir(fullfile(param.SaveDirectory))
    end
            
                 
    Utr = param.WB.Utr;
    S = param.WB.S;

    save(fullfile(param.SaveDirectory,['WB_Utr_',num2str(param.constW),'.mat']), 'Utr','-v7.3')
    save(fullfile(param.SaveDirectory,['WB_S_',num2str(param.constW),'.mat']), 'S','-v7.3')
 
end