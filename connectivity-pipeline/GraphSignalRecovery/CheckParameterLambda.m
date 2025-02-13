
function lambda  = CheckParameterLambda(param,BB,L,BC)

    lambdacheck = fullfile(param.SaveInpainting, ['RegularizationParameters_LCurve_',param.subject,'.mat']);
    
    %if ~exist(lambdacheck, 'file')
        
        % Load preliminary fMRI volume
        Inputs_fMRI
        
        i = 20; % frame to evaluate
        vol = V(i,:);
        
        clear V
        
        disp('Calculating correct parameter lambda..')
        lambdas = 1:100;
%         diffs = zeros(length(lambdas),1);
        x1norms = zeros(length(lambdas),1);
        for c = 1:length(lambdas)
            lambda = lambdas(c);
            L_mat = BB+lambda*L;
            L1 = ichol(L_mat);
%             alpha = max(sum(abs(L_mat),2)./diag(L_mat))-2;
%             L1 = ichol(L_mat, struct('type','ict','droptol',1e-3,'diagcomp',alpha));
            xsignal = zeros(size(L1,1),1);
            xsignal(BC) = vol(BC);
            [x1,~,~,~,~] = pcg(L_mat,xsignal,1e-8,100,L1,L1');
            x1norms(c) = norm(x1'*(L*x1));
%             diff = B*x1-xsignal;
%             diffs(c) = norm(diff)^2;
        end

        lambdaStars.lambdas = lambdas;
        lambdaStars.x1norms = x1norms;
        save(fullfile(param.SaveInpainting, ['RegularizationParameters_LCurve_',param.subject,'.mat']),'lambdaStars')
        
%     else
%         load(lambdacheck)
%         x1norms = lambdaStars.x1norms;
%         lambdas = lambdaStars.lambdas;
    %end
    
    lambda = knee_pt(x1norms,lambdas);
end