clear
addpath GA/
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Mutation= getAvailableMutation(prepro_BRCA);
Mutation= rmirrelevant(Mutation);

K = 5;
N = length(Mutation.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
    for i = 1:K
        Train.F=Mutation.Features(:,Folds ~= i); 
        Train.S=Mutation.Survival(Folds ~= i); 
        Train.C=Mutation.Censored(Folds ~= i); 
        Test.F=Mutation.Features(:,Folds == i); 
        Test.S=Mutation.Survival(Folds == i); 
        Test.C=Mutation.Censored(Folds == i); 
   
        p = size(Train.F,1); % p=number of feature
        GenomeLength = p; % This is the number of features in the dataset
        options = gaoptimset('CreationFcn', {@PopFunction},...
                             'PopulationSize',10,...
                             'Generations',10,...
                             'PopulationType', 'bitstring',... 
                             'SelectionFcn',{@selectiontournament,2},...
                             'MutationFcn',{@mutationuniform, 0.1},...
                             'CrossoverFcn', {@crossoverarithmetic,0.8},...
                             'EliteCount',2,...
                             'StallGenLimit',10,...
                             'PlotFcns',{@gaplotbestf},...  
                             'Display', 'iter',...
                             'UseParallel', true); 
          
        FitnessFunction = @(x)c_index_fitness(x, Train);
        [chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
        Best_feature_set = chromosome;% Best feature set
        Best_feature_Index = find(Best_feature_set==1) % Index of Chromosome

        Beta = coxphfit(Train.F(Best_feature_Index,:).',  Train.S.','Censoring',  Train.C.');
        C(i) = cIndex(Beta,  Test.F(Best_feature_Index,:).', Test.S, Test.C);
        
    end
    
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds C Basic Beta i K N original_BRCA prepro_BRCA