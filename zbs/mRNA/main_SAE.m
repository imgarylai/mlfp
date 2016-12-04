clear
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
mRNA= getAvailablemRNA(prepro_BRCA);
mRNA= rmirrelevant(mRNA);
addpath sparse_autoencoder/
addpath sparse_autoencoder/minFunc
% Report testing error using CI and 5 trial 
K = 5;
N = length(mRNA.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
for i = 1:K
    Basic=mRNA.Features;
    %%%%%%SAE%%%%%%
    [ TrainFeature, mRNA_weight,b ] = sparse_autoencoder(Basic(:, Folds ~= i), 10);
    TestFeature = s_encode(mRNA_weight,b, Basic(:, Folds == i));
    %%%%%%%%%%%%%%%
    Beta = coxphfit(TrainFeature.', mRNA.Survival(Folds ~= i).',...
        'Censoring', mRNA.Censored(Folds ~= i).');
    C(i) = cIndex(Beta, TestFeature.', mRNA.Survival(Folds == i),...
        mRNA.Censored(Folds == i));
    clear TrainFeature TestFeature
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA