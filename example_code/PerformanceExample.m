%Demonstration of performance assessment techniques using class project
%dataset.

% turn off warnings
warning('off','all')

% read in data and extract minimal features
load '../data/GBMLGG.Data.mat';

%define a basic model containing IDH mutations, age and chromosomes 1p/19q
x1 = Features(strcmp(cellstr(Symbols),...
    'age_at_initial_pathologic_diagnosis_Clinical'), :);
x2 = Features(strcmp(cellstr(Symbols), 'IDH1_Mut'), :);
x3 = Features(strcmp(cellstr(Symbols), 'IDH2_Mut'), :);
x4 = Features(strcmp(cellstr(Symbols), '1p_CNVArm'), :);
x5 = Features(strcmp(cellstr(Symbols), '19q_CNVArm'), :);
Basic = [x1; x2; x3; x4; x5];

%convert symbols, symbol types to cell array
Symbols = cellstr(Symbols);
SymbolTypes = cellstr(SymbolTypes);

%define core set of samples that have all basic features, mutation features,
%survival, censoring
Keep = ~isnan(Survival) & ~isnan(Censored) ...
    & (sum(isnan(Basic), 1) == 0) ...
    & (sum(isnan(Features(strcmp(SymbolTypes, 'CNVArm'),:)), 1) == 0);
Basic = Basic(:, Keep);
Features = Features(strcmp(SymbolTypes, 'CNVArm'), Keep);
Symbols = Symbols(strcmp(SymbolTypes, 'CNVArm'));
Survival = Survival(Keep);
Censored = Censored(Keep);
N = length(Survival);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard train/test/validation - model selection and assessment %%%%%%%%%

fprintf('Model assessment and selection with train/test/validation approach.\n');

%define train, test, and validation sets
Train = [1 : ceil(N/2)];
Validation = [ceil(N/2)+1 : ceil(3*N/4)];
Test = [ceil(3*N/4)+1 : N];

%define performance of basic model
BetaBasic = coxphfit(Basic(:, Train).', Survival(Train).',...
    'Censoring', Censored(Train).');
cBasic = cIndex(BetaBasic, Basic(:, Validation).', Survival(Validation),...
    Censored(Validation));

%perform 1 iteration of forward stagewise feature selection to grow model
Improvement = zeros(size(Features, 1), 1);
for i = 1:size(Features, 1)
    
    %extract candidate feature values
    NewTrain = Features(i, Train);
    NewVal = Features(i, Validation);
    
    %if feature is not too sparse then build a model and test improvement
    if ~isnan(var(NewTrain))
        Beta = coxphfit([Basic(:, Train); NewTrain].',...
            Survival(Train).', 'Censoring', Censored(Train).');
        cAugmented = cIndex(Beta, [Basic(:, Validation); NewVal].',...
            Survival(Validation), Censored(Validation));
        Improvement(i) = cAugmented - cBasic;
    end
    
end
[~, Best] = max(Improvement);

%retrain model on train+validation samples and test
BetaStar = coxphfit([Basic(:, [Train Validation]); Features(Best, [Train Validation])].',...
    Survival([Train Validation]).', 'Censoring', Censored([Train Validation]).');
C = cIndex(BetaStar, [Basic(:, Test);  Features(Best, Test)].', Survival(Test), Censored(Test));

%report performance of best model
fprintf('\tModel - with %s: c-index = %g\n', Symbols{Best}, C);

%Build final model with complete data - this is the model that you deploy.
%The predicted generalization error of this model is 'C'
BetaFinal = coxphfit([Basic; Features(Best, :)].',...
    Survival.', 'Censoring', Censored.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10-fold cross validation of a basic model (no model selection) %%%%%%%%%%

%define number of folds and assign samples
K = 10;
Folds = ceil([1:N] / (N/K));

fprintf('Model assessment with %d-fold cross-validation.\n', K);

%cycle through folds, training and testing model
C = nan(1,K);
for i = 1:K
    Beta = coxphfit(Basic(:, Folds ~= i).', Survival(Folds ~= i).',...
        'Censoring', Censored(Folds ~= i).');
    C(i) = cIndex(Beta, Basic(:, Folds == i).', Survival(Folds == i),...
        Censored(Folds == i));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested 5-fold cross validation - model selection and assessment %%%%%%%%%

fprintf('Nested cross validation.\n');

%define train + validation (75%), test (25%) sets
Train = [1 : ceil(3*N/4)];
Test = [ceil(3*N/4)+1 : N];

%define number of folds and assign samples in training + validation
K = 5;
Folds = ceil([1:ceil(3*N/4)] / (ceil(3*N/4)/K));

%perform 1-iteration of forward stagewise feature selection to improve model
cAugmented = nan(size(Features, 1), K);
for i = 1:size(Features, 1)
    for j = 1:K
        
        %extract candidate feature values
        NewTrain = Features(i, Train(Folds ~= j));
        NewVal = Features(i, Train(Folds == j));
        
        %if feature is not too sparse then build a model and test improvement
        if ~isnan(var(NewTrain))
            Beta = coxphfit([Basic(:, Train(Folds ~= j)); NewTrain].',...
                Survival(Train(Folds ~= j)).',...
                'Censoring', Censored(Train(Folds ~= j)).');
            cAugmented(i,j) = cIndex(Beta, [Basic(:, Train(Folds == j)); NewVal].',...
                Survival(Train(Folds == j)), Censored(Train(Folds == j)));
        end
    end
end

%pick the feature with the best average improvement
Average = mean(cAugmented, 2);
[Max, Best] = max(Average);

%retrain model on all data and test
BetaStar = coxphfit([Basic(:, Train); Features(Best, Train)].',...
    Survival(Train).', 'Censoring', Censored(Train).');
C = cIndex(BetaStar, [Basic(:, Test); Features(Best, Test)].',...
            Survival(Test), Censored(Test));

%report performance of best model on testing data
fprintf('\tNested CV model - with %s: c-index = %g\n', Symbols{Best}, C);

%Build final model with complete data - this is the model that you deploy.
%The predicted generalization error of this model is 'C'
BetaFinal = coxphfit([Basic; Features(Best, :)].',...
    Survival.', 'Censoring', Censored.');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Nested 10-fold cross validation - repeated trials %%%%%%%%%%%%%%%%%%%%%%%

% define number of trials and initialize testing error container
Trials = 10;
C = zeros(Trials, 1); %testing c-index
Best = zeros(Trials, 1); %index of best feature

fprintf('Repeated trials nested cross validation.\n');

for t = 1:Trials
    
    %generate a shuffling of the samples
    [~, Shuffle] = sort(rand(N, 1));
    
    %define train + validation (75%), test (25%) sets
    Train = Shuffle(1 : ceil(3*N/4));
    Test = Shuffle(ceil(3*N/4)+1 : N);

    %define number of folds and assign samples in training + validation
    K = 5;
    Folds = ceil([1:ceil(3*N/4)] / (ceil(3*N/4)/K));
    
    %perform 1-iteration of forward stagewise feature selection to improve model
    Improvement = zeros(size(Features, 1), K);
    for i = 1:size(Features, 1)
        
        %extract candidate feature values
        NewTrain = Features(i, Train(Folds ~= j));
        NewVal = Features(i, Train(Folds == j));
        
        %if feature is not too sparse then build a model and test improvement
        if ~isnan(var(NewTrain))
            Beta = coxphfit([Basic(:, Train(Folds ~= j)); NewTrain].',...
                Survival(Train(Folds ~= j)).',...
                'Censoring', Censored(Train(Folds ~= j)).');
            cAugmented(i,j) = cIndex(Beta, [Basic(:, Train(Folds == j)); NewVal].',...
                Survival(Train(Folds == j)), Censored(Train(Folds == j)));
        end
    end
    
    %pick the feature with the best average improvement
    Average = mean(cAugmented, 2);
    [Max, Best(t)] = max(Average);
    
    %retrain model on all data and test
    BetaStar = coxphfit([Basic(:, Train); Features(Best(t), Train)].',...
                        Survival(Train).', 'Censoring', Censored(Train).');
    C(t) = cIndex(BetaStar, [Basic(:, Test); Features(Best(t), Test)] .',...
                     Survival(Test), Censored(Test));
    
end

%report performance of best model on testing data
fprintf('\tNested CV model with %d shuffles: average c-index = %g\n',...
        Trials, mean(C));

%Build final model with complete data - here you can take a couple
%approaches:
%1 - Use the feature that provided the best error. This is not cheating 
%    because you are reporting the average error - not the error
%    corresponding to this model.
%2 - Use the feature that was selected most often in the experiments. If
%    there is not a clear winner then your data is not stable and is very
%    sensitive to partitioning into training / testing.

%approach 1
[~, Optimal] = max(C);
BetaFinal1 = coxphfit([Basic; Features(Best(Optimal), :)].',...
    Survival.', 'Censoring', Censored.');

%approach 2
Optimal = mode(Best);
BetaFinal2 = coxphfit([Basic; Features(Optimal, :)].',...
    Survival.', 'Censoring', Censored.');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .632 Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Model assessment with 0.632 boostrap.\n');

%number of bootstrap samples
B = 100;

%build model on full dataset, calculate training c-index
Beta = coxphfit(Basic.', Survival, 'Censoring', Censored );
err = cIndex(Beta, Basic.', Survival, Censored);

%generate boostrap samples and track error
Cboot = nan(1,B);
for i = 1:B
    
    %generate random sample index
    Sample = ceil(length(Survival) * rand(1,length(Survival)));
    while(sum(Censored(Sample)) == length(Sample))
        Sample = ceil(length(Survival) * rand(1,length(Survival)));
    end
    
    %generate out of bootstrap sample
    Out = setdiff([1:N], Sample);
    
    %build model
    Beta = coxphfit(Basic(:, Sample).', Survival(Sample),...
        'censoring', Censored(Sample));
    
    %compute 'c' for bootstrap model on left out samples
    Cboot(i) = cIndex(Beta, Basic(:, Out).', Survival(Out), Censored(Out));
    
end
C = 0.368 * err + 0.632 * sum(Cboot) / B;
fprintf('\tTraining = %g, Bootstrap = %g, 0.632 c-index = %g\n',...
    err, sum(Cboot)/B, C);


% turn warnings back on
warning('on','all')