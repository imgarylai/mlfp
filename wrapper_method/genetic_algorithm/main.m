clear
warning('off', 'all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA = preprocessing(original_BRCA);
[train,test] = getTrainigAndTesting(prepro_BRCA);
clinical = getAvailableClinical(train);
mutation = getAvailableMutation(train);
CNV = getAvailableCNV(train);
protein = getAvailableProtein(train);

p = size(clinical.Features,1); % p=number of feature
GenomeLength = p; % This is the number of features in the dataset

options = gaoptimset('CreationFcn', {@PopFunction},...
                     'PopulationSize',100,...
                     'Generations',1000,...
                     'PopulationType', 'bitstring',...
                     'SelectionFcn',{@selectiontournament,4},...
                     'MutationFcn',{@mutationuniform, 0.1},...
                     'CrossoverFcn', {@crossoverarithmetic,0.8},...
                     'EliteCount',2,...
                     'StallGenLimit',200,...
                     'PlotFcns',{@gaplotbestf},...
                     'Display', 'iter',...
                     'UseParallel', true);

%'UseParallel', true
nVars = p;
FitnessFunction = @(x)c_index_fitness(x, clinical);

[chromosome, ~, ~, ~, ~, ~] = ga(FitnessFunction, nVars, options);
Best_feature_set = chromosome % Best feature set
Best_feature_Index = find(Best_feature_set == 1) % Index of Chromosome
