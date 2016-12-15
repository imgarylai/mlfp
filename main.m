clear
addpath autoencoder/
addpath genetic_algorithm/
addpath sparse_autoencoder/
addpath sparse_autoencoder/minFunc/
addpath PCA/
warning('off','all')
prompt_dateset = '1->BRCA\n2->GBMLGG\n (enter a number and press enter: )\n';
dateset = input(prompt_dateset);
prompt_datatype = '1->Clinical\n2->Mutation\n3->CNV\n4->Protein\n5->mRNA\n(enter a number and press enter: )\n';
datatype = input(prompt_datatype);
if (datatype==1 | datatype==2 | datatype==3|datatype==4)
    prompt_method = '1->no\n2->GA\n3->AE\n4->S-AE\n5->PCA\n6->AE+GA\n7->S-AE+GA\n8->PCA+GA\n(enter a number and press enter: )\n';
    method = input(prompt_method);
elseif (datatype==5)  
    prompt_method = '3->AE\n4->S-AE\n5->PCA\n6->AE+GA\n7->S-AE+GA\n8->PCA+GA\n(enter a number and press enter: )\n';
    method = input(prompt_method);
end
if (method==3|method==4|method==5|method==6|method==7|method==8)
    prompt_dimension = 'enter the dimension:\n';
    dimension = input(prompt_dimension);
end

if (dateset==1)
    load ('data/BRCA.mat');
    original = BRCA;
elseif (dateset==2)
    load ('data/GBMLGG.mat');
    original = GBMLGG;
end

if (datatype==1)
    Data= original.Clinical;
elseif (datatype==2)
    Data=original.Mutation;
elseif (datatype==3)
    Data= original.CNV;
elseif (datatype==4)
    Data= original.Protein;
elseif (datatype==5)    
    Data= original.mRNA;   
end

 K = 5;
 N = length(Data.Survival);
 Folds = ceil([1:N] / (N/K));
 C = nan(1,K);


if (method==1)   % No feature selection, no dimension deduction
    for i = 1:K
        Basic=Data.Features; 
        TrainFeature = Basic(:, Folds ~= i);
        TestFeature  = Basic(:, Folds == i);
        Beta = coxphfit(TrainFeature.', Data.Survival(Folds ~= i).',...
            'Censoring', Data.Censored(Folds ~= i).');
        C(i) = cIndex(Beta,  TestFeature .', Data.Survival(Folds == i),...
            Data.Censored(Folds == i));
        clear TrainFeature TestFeature
    end
end  

if (method==2)  % Genetic Algorithm 
    for i = 1:K
        Train.F=Data.Features(:,Folds ~= i); 
        Train.S=Data.Survival(Folds ~= i); 
        Train.C=Data.Censored(Folds ~= i); 
        Test.F=Data.Features(:,Folds == i); 
        Test.S=Data.Survival(Folds == i); 
        Test.C=Data.Censored(Folds == i); 
   
        p = size(Train.F,1); % p=number of feature
        GenomeLength = p; % This is the number of features in the dataset
        options = gaoptimset('CreationFcn', {@PopFunction},...
                             'PopulationSize',30,...
                             'Generations',50,...
                             'PopulationType', 'bitstring',... 
                             'SelectionFcn',{@selectiontournament,2},...
                             'MutationFcn',{@mutationuniform, 0.1},...
                             'CrossoverFcn', {@crossoverarithmetic,0.8},...
                             'EliteCount',2,...
                             'StallGenLimit',100,...
                             'PlotFcns',{@gaplotbestf},...  
                             'Display', 'iter',...
                             'UseParallel', true); 
          
        FitnessFunction = @(x)c_index_fitness(x, Train);
        [chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
        Best_feature_set = chromosome;% Best feature set
        Final_feature_set(i,:) = Best_feature_set;       
    end
        Best=sum(Final_feature_set);
        Best_feature_Index=find(Best>=4)
        ga_dimension=length(Best_feature_Index);
        for j = 1:K
        Train.F=Data.Features(:,Folds ~= j); 
        Train.S=Data.Survival(Folds ~= j); 
        Train.C=Data.Censored(Folds ~= j); 
        Test.F=Data.Features(:,Folds == j); 
        Test.S=Data.Survival(Folds == j); 
        Test.C=Data.Censored(Folds == j); 
        Beta = coxphfit(Train.F(Best_feature_Index,:).',  Train.S.','Censoring',  Train.C.');
        C(j) = cIndex(Beta,  Test.F(Best_feature_Index,:).', Test.S, Test.C);      
        end
end  

if (method==3) % Autoencoder
        for i = 1:K
            Basic=Data.Features; 
            [TrainFeature, Data_weight] = autoencoder(Basic(:, Folds ~= i), dimension);
            TestFeature = encode(Data_weight, Basic(:, Folds == i));
            Beta = coxphfit(TrainFeature.', Data.Survival(Folds ~= i).',...
                'Censoring', Data.Censored(Folds ~= i).');
            C(i) = cIndex(Beta,  TestFeature .', Data.Survival(Folds == i),...
                Data.Censored(Folds == i));
            clear TrainFeature TestFeature Beta
        end
end  


if (method==4) 
    % Sparse-Autoencoder
    for i = 1:K
        Basic=Data.Features;
        [ TrainFeature, Data_weight,b ] = sparse_autoencoder(Basic(:, Folds ~= i), dimension);
        TestFeature = s_encode(Data_weight,b, Basic(:, Folds == i));  
        Beta = coxphfit(TrainFeature.', Data.Survival(Folds ~= i).',...
            'Censoring', Data.Censored(Folds ~= i).');
        C(i) = cIndex(Beta, TestFeature.', Data.Survival(Folds == i),...
            Data.Censored(Folds == i));
        clear TrainFeature TestFeature
    end
    
end  

if (method==5)   % PCA
    for i = 1:K
        Basic=Data.Features;

        [TrainFeature, mapping] = newPCA(Basic(:, Folds ~= i)', dimension);        
        TestFeature= Basic(:, Folds == i)' * mapping.M; 
        TestFeature = TestFeature - mean(TestFeature, 1); 
     
        TrainFeature=TrainFeature';
        TestFeature=TestFeature';

        Beta = coxphfit(TrainFeature.', Data.Survival(Folds ~= i).',...
            'Censoring', Data.Censored(Folds ~= i).');
        C(i) = cIndex(Beta, TestFeature.', Data.Survival(Folds == i),...
            Data.Censored(Folds == i));
        clear TrainFeature TestFeature
    end
end  


if (method==6)   % AE+GA
    for i = 1:K
        Basic=Data.Features; 
        [TrainFeature, Data_weight] = autoencoder(Basic(:, Folds ~= i), dimension);
        TestFeature = encode(Data_weight, Basic(:, Folds == i));
        eval(['TrainFeature' num2str(i) '= TrainFeature']);
        eval(['TestFeature' num2str(i) '= TestFeature']);
        
        Train.F=TrainFeature; 
        Train.S=Data.Survival(Folds ~= i); 
        Train.C=Data.Censored(Folds ~= i); 
        Test.F=TestFeature; 
        Test.S=Data.Survival(Folds == i); 
        Test.C=Data.Censored(Folds == i); 
        p = size(Train.F,1); % p=number of feature
        GenomeLength = p; % This is the number of features in the dataset
        options = gaoptimset('CreationFcn', {@PopFunction},...
                             'PopulationSize',30,...
                             'Generations',50,...
                             'PopulationType', 'bitstring',... 
                             'SelectionFcn',{@selectiontournament,2},...
                             'MutationFcn',{@mutationuniform, 0.1},...
                             'CrossoverFcn', {@crossoverarithmetic,0.8},...
                             'EliteCount',2,...
                             'StallGenLimit',100,...
                             'PlotFcns',{@gaplotbestf},...  
                             'Display', 'iter',...
                             'UseParallel', true); 
          
        FitnessFunction = @(x)c_index_fitness(x, Train);
        [chromosome,~,~,~,~,~] = ga(FitnessFunction,p,options);
        Best_feature_set = chromosome;% Best feature set
        Final_feature_set(i,:) = Best_feature_set;       
    end
    Best=sum(Final_feature_set);
    Best_feature_Index=find(Best>=4)
    ga_dimension=length(Best_feature_Index);
        for j = 1:K
%         Basic=Data.Features; 
%         [TrainFeature, Data_weight] = autoencoder(Basic(:, Folds ~= j), dimension);
%         TestFeature = encode(Data_weight, Basic(:, Folds == j));
        Train.F=eval(['TrainFeature' num2str(j)]); 
        Train.S=Data.Survival(Folds ~= j); 
        Train.C=Data.Censored(Folds ~= j); 
        Test.F=eval(['TestFeature' num2str(j)]);  
        Test.S=Data.Survival(Folds == j); 
        Test.C=Data.Censored(Folds == j); 
        Beta = coxphfit(Train.F(Best_feature_Index,:).',  Train.S.','Censoring',  Train.C.');
        C(j) = cIndex(Beta,  Test.F(Best_feature_Index,:).', Test.S, Test.C);      
        end 
end  













if(method ==1)
    fprintf('\tDimension: %g --> %g \n', size(Data.Features,1), size(Data.Features,1));
elseif (method ==2)
    fprintf('\tDimension: %g --> %g \n', size(Data.Features,1), mean(ga_dimension));
elseif (method==3|method==4|method==5)
    fprintf('\tDimension: %g --> %g\n', size(Data.Features,1), dimension);
elseif (method==6|method==7|method==8)
    fprintf('\tDimension: %g --> %g --> %g\n', size(Data.Features,1), dimension,mean(ga_dimension));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));
