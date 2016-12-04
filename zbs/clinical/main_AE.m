clear
addpath autoencoder/
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Clinical= getAvailableClinical(prepro_BRCA);
Clinical= rmirrelevant(Clinical);

K = 5;
N = length(Clinical.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
    for i = 1:K
        Basic=Clinical.Features; 
        [TrainFeature, Clinical_weight] = autoencoder(Basic(:, Folds ~= i), 20);
         TestFeature = encode(Clinical_weight, Basic(:, Folds == i));
%         TrainFeature = Basic(:, Folds ~= i);
%         TestFeature  = Basic(:, Folds == i);

        Beta = coxphfit(TrainFeature.', Clinical.Survival(Folds ~= i).',...
            'Censoring', Clinical.Censored(Folds ~= i).');
        C(i) = cIndex(Beta,  TestFeature .', Clinical.Survival(Folds == i),...
            Clinical.Censored(Folds == i));
        clear TrainFeature TestFeature
    end
    
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA