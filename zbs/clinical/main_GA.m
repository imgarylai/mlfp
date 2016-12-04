clear
addpath GA/
warning('off','all')
original_BRCA = load ('BRCA.Data.mat');
prepro_BRCA=prepro(original_BRCA);
Clinical= getAvailableClinical(prepro_BRCA);
Clinical= rmirrelevant(Clinical);
K = 5;
N = length(Clinical.Survival);
Folds = ceil([1:N] / (N/K));
C = nan(1,K);
    for i = 1:1
        Train.F=Clinical.Features(:,Folds ~= i); 
        Train.S=Clinical.Survival(Folds ~= i); 
        Train.C=Clinical.Censored(Folds ~= i); 
        Test.F=Clinical.Features(:,Folds == i); 
        Test.S=Clinical.Survival(Folds == i); 
        Test.C=Clinical.Censored(Folds == i); 
       
        p = size(Train.F,1); % p=number of feature
        GenomeLength = p; % This is the number of features in the dataset
        options = gaoptimset('CreationFcn', {@PopFunction},...
                             'PopulationSize',20,...
                             'Generations',300,...
                             'PopulationType', 'bitstring',... 
                             'SelectionFcn',{@selectiontournament,2},...
                             'MutationFcn',{@mutationuniform, 0.1},...
                             'CrossoverFcn', {@crossoverarithmetic,0.8},...
                             'EliteCount',2,...
                             'StallGenLimit',300,...
                             'PlotFcns',{@gaplotbestf},...  
                             'Display', 'iter',...
                             'UseParallel', true); 
          
        FitnessFunction = @(x)c_index_fitness(x,Train);
        [chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
        Best_feature_set = chromosome; % Best feature set
        Best_feature_Index = find(Best_feature_set==1) % Index of Chromosome
        Beta = coxphfit(Train.F(Best_feature_Index,:).',  Train.S.','Censoring',  Train.C.');
        C(i) = cIndex(Beta,  Test.F(Best_feature_Index,:).', Test.S, Test.C);
    end
 
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
clear Folds Basic Beta i K N original_BRCA prepro_BRCA