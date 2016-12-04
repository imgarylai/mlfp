clear
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Mutation= getAvailableMutation(prepro_BRCA);
Mutation= rmirrelevant(Mutation);
    addpath GA/

% Report testing error using CI and 5 trial 
K = 5;
N = length(Mutation.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
for i = 1:1
    Basic=Mutation.Features;
    [ Mutation.TrainFeature, Mutation_weight ] = autoencoder(Basic(:, Folds ~= i), 30);
      Mutation.TestFeature = encode(Mutation_weight, Basic(:, Folds == i));
     Mutation.trainSurvival=Mutation.Survival(Folds ~= i);
     Mutation.trainCensored=Mutation.Censored(Folds ~= i);

p = size(Mutation.TrainFeature,1); % p=number of feature
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
          
FitnessFunction = @(x)c_index_fitness(x, Mutation);
[chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
Best_feature_set = chromosome;% Best feature set
Best_feature_Index = find(Best_feature_set==1) % Index of Chromosome


    Beta = coxphfit(Mutation.TrainFeature([Best_feature_Index],:).', Mutation.trainSurvival(:).',...
        'Censoring', Mutation.trainCensored(:).');
    C(i) = cIndex(Beta, Mutation.TestFeature([Best_feature_Index],:).', Mutation.Survival(Folds == i),...
        Mutation.Censored(Folds == i));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA






