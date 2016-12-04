clear
addpath autoencoder/

warning('off','all')
original_BRCA = load ('/Users/Bashan/Desktop/mymlfp/data/BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Mutation= getAvailableMutation(prepro_BRCA);
Mutation= rmirrelevant(Mutation);
% 
% K = 5;
% N = length(Mutation.Survival);
% Folds = ceil([1:N] / (N/K));
% C = nan(1,K);
%     for i = 1:K
%         Basic=Mutation.Features; 
%         [TrainFeature, Mutation_weight] = autoencoder(Basic(:, Folds ~= i), 20);
%          TestFeature = encode(Mutation_weight, Basic(:, Folds == i));
% %         TrainFeature = Basic(:, Folds ~= i);
% %         TestFeature  = Basic(:, Folds == i);
% 
%         Beta = coxphfit(TrainFeature.', Mutation.Survival(Folds ~= i).',...
%             'Censoring', Mutation.Censored(Folds ~= i).');
%         C(i) = cIndex(Beta,  TestFeature .', Mutation.Survival(Folds == i),...
%             Mutation.Censored(Folds == i));
%         clear TrainFeature TestFeature
%     end
%     
% fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
% clear Folds C Basic Beta i K N original_BRCA prepro_BRCA