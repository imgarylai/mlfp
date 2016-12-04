clear
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
CNV= getAvailableCNV(prepro_BRCA);
CNV= rmirrelevant(CNV);
    addpath GA/

% Report testing error using CI and 5 trial 
K = 5;
N = length(CNV.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
for i = 1:1
    Basic=CNV.Features;
    [ CNV.TrainFeature, CNV_weight ] = autoencoder(Basic(:, Folds ~= i), 30);
      CNV.TestFeature = encode(CNV_weight, Basic(:, Folds == i));
     CNV.trainSurvival=CNV.Survival(Folds ~= i);
     CNV.trainCensored=CNV.Censored(Folds ~= i);

p = size(CNV.TrainFeature,1); % p=number of feature
GenomeLength = p; % This is the number of features in the dataset
options = gaoptimset('CreationFcn', {@PopFunction},...
                     'PopulationSize',30,...
                     'Generations',100,...
                     'PopulationType', 'bitstring',... 
                     'SelectionFcn',{@selectiontournament,2},...
                     'MutationFcn',{@mutationuniform, 0.2},...
                     'CrossoverFcn', {@crossoverarithmetic,0.7},...
                     'EliteCount',2,...
                     'StallGenLimit',100,...
                     'PlotFcns',{@gaplotbestf},...  
                     'Display', 'iter',...
                     'UseParallel', true); 
          
FitnessFunction = @(x)c_index_fitness(x, CNV);
[chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
Best_feature_set = chromosome;% Best feature set
Best_feature_Index = find(Best_feature_set==1) % Index of Chromosome


    Beta = coxphfit(CNV.TrainFeature([Best_feature_Index],:).', CNV.trainSurvival(:).',...
        'Censoring', CNV.trainCensored(:).');
    C(i) = cIndex(Beta, CNV.TestFeature([Best_feature_Index],:).', CNV.Survival(Folds == i),...
        CNV.Censored(Folds == i));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA






