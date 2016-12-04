clear
addpath autoencoder/

warning('off','all')
original_BRCA = load ('/Users/Bashan/Desktop/mymlfp/data/BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
CNV= getAvailableCNV(prepro_BRCA);
CNV= rmirrelevant(CNV);
% 
% K = 5;
% N = length(CNV.Survival);
% Folds = ceil([1:N] / (N/K));
% C = nan(1,K);
%     for i = 1:K
%         Basic=CNV.Features; 
%         [TrainFeature, CNV_weight] = autoencoder(Basic(:, Folds ~= i), 20);
%          TestFeature = encode(CNV_weight, Basic(:, Folds == i));
% %         TrainFeature = Basic(:, Folds ~= i);
% %         TestFeature  = Basic(:, Folds == i);
% 
%         Beta = coxphfit(TrainFeature.', CNV.Survival(Folds ~= i).',...
%             'Censoring', CNV.Censored(Folds ~= i).');
%         C(i) = cIndex(Beta,  TestFeature .', CNV.Survival(Folds == i),...
%             CNV.Censored(Folds == i));
%         clear TrainFeature TestFeature
%     end
%     
% fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
% clear Folds C Basic Beta i K N original_BRCA prepro_BRCA