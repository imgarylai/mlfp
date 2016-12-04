clear
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Protein= getAvailableProtein(prepro_BRCA);
Protein= rmirrelevant(Protein);
    addpath GA/

% Report testing error using CI and 5 trial 
K = 5;
N = length(Protein.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
for i = 1:1
    Basic=Protein.Features;
    [ Protein.TrainFeature, Protein_weight ] = autoencoder(Basic(:, Folds ~= i), 30);
      Protein.TestFeature = encode(Protein_weight, Basic(:, Folds == i));
     Protein.trainSurvival=Protein.Survival(Folds ~= i);
     Protein.trainCensored=Protein.Censored(Folds ~= i);

p = size(Protein.TrainFeature,1); % p=number of feature
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
          
FitnessFunction = @(x)c_index_fitness(x, Protein);
[chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
Best_feature_set = chromosome;% Best feature set
Best_feature_Index = find(Best_feature_set==1) % Index of Chromosome


    Beta = coxphfit(Protein.TrainFeature([Best_feature_Index],:).', Protein.trainSurvival(:).',...
        'Censoring', Protein.trainCensored(:).');
    C(i) = cIndex(Beta, Protein.TestFeature([Best_feature_Index],:).', Protein.Survival(Folds == i),...
        Protein.Censored(Folds == i));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA






