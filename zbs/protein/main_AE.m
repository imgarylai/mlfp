clear
addpath autoencoder/

warning('off','all')
original_BRCA = load ('/Users/Bashan/Desktop/mymlfp/data/BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Protein= getAvailableProtein(prepro_BRCA);
Protein= rmirrelevant(Protein);
% 
% K = 5;
% N = length(Protein.Survival);
% Folds = ceil([1:N] / (N/K));
% C = nan(1,K);
%     for i = 1:K
%         Basic=Protein.Features; 
%         [TrainFeature, Protein_weight] = autoencoder(Basic(:, Folds ~= i), 20);
%          TestFeature = encode(Protein_weight, Basic(:, Folds == i));
% %         TrainFeature = Basic(:, Folds ~= i);
% %         TestFeature  = Basic(:, Folds == i);
% 
%         Beta = coxphfit(TrainFeature.', Protein.Survival(Folds ~= i).',...
%             'Censoring', Protein.Censored(Folds ~= i).');
%         C(i) = cIndex(Beta,  TestFeature .', Protein.Survival(Folds == i),...
%             Protein.Censored(Folds == i));
%         clear TrainFeature TestFeature
%     end
%     
% fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
% clear Folds C Basic Beta i K N original_BRCA prepro_BRCA