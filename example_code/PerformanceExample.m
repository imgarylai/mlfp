%Demonstration of performance assessment techniques using class project
%dataset.

% read in data and extract minimal features
load './GBMLGG.Data.mat';
x1 = Features(strcmp(cellstr(Symbols),...
                       'age_at_initial_pathologic_diagnosis_Clinical'), :);
x2 = Features(strcmp(cellstr(Symbols), 'IDH1_Mut'), :);
x3 = Features(strcmp(cellstr(Symbols), 'IDH2_Mut'), :);
x4 = Features(strcmp(cellstr(Symbols), '1p_CNVArm'), :);
x5 = Features(strcmp(cellstr(Symbols), '19q_CNVArm'), :);
X = [x1; x2; x3; x4; x5];

% define core set of samples that have all features, survival, censoring
Keep = ~isnan(Survival) & ~isnan(Censored) & (sum(isnan(X), 1) == 0);
X = X(:, Keep);
Survival = Survival(Keep);
Censored = Censored(Keep);
N = length(Survival);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard train/test/validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Model assessment and selection with train/test/validation approach.\n');

%simulate feature selection by comparing two models: with/without NOTCH1
NOTCH1 = Features(strcmp(cellstr(Symbols), 'NOTCH1_Mut'), Keep);

%define train, test, and validation sets
Train = [1 : ceil(N/2)];
Validation = [ceil(N/2)+1 : ceil(3*N/4)];
Test = [ceil(3*N/4)+1 : N];

%build each model on training
Beta1 = coxphfit(X(:, Train).', Survival(Train).',...
                 'Censoring', Censored(Train).');
Beta2 = coxphfit([X(:, Train); NOTCH1(Train)].', Survival(Train).',...
                 'Censoring', Censored(Train).');

%compare performance of models on validation
C1 = cIndex(Beta1, X(:, Validation).', Survival(Validation),...
    Censored(Validation));
C2 = cIndex(Beta2, [X(:, Validation); NOTCH1(Validation)].',...
    Survival(Validation), Censored(Validation));

%report performance of best model on validation
if(C1 >= C2)
    C = cIndex(Beta1, X(:, Test).', Survival(Test), Censored(Test));
    fprintf('\tModel 1 - without NOTCH1: c-index = %g\n', C);
else
    C = cIndex(Beta2, [X(:, Test); NOTCH1(Test)].', Survival(Test),...
        Censored(Test));
    fprintf('\tModel 2 - with NOTCH1: c-index = %g\n', C);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 10-fold cross validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define number of folds and assign samples
K = 10;
Folds = ceil([1:N]/(N/K));

fprintf('Model assessment with %d-fold cross-validation.\n', K);

%cycle through folds, training and testing model
C = nan(1,K);
for i = 1:K
    Beta = coxphfit(X(:, Folds ~= i).', Survival(Folds ~= i).',...
                 'Censoring', Censored(Folds ~= i).');
    C(i) = cIndex(Beta, X(:, Folds == i).', Survival(Folds == i),...
                 Censored(Folds == i));
end
fprintf('\tmean c-index = %g, standard deviation = %g\n', mean(C), std(C));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% .632 Bootstrap %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Model assessment with 0.632 boostrap.\n');

%number of bootstrap samples
B = 100;

%build model on full dataset, calculate training c-index
Beta = coxphfit(X.', Survival, 'Censoring', Censored );
err = cIndex(Beta, X.', Survival, Censored);

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
    Beta = coxphfit(X(:, Sample).', Survival(Sample),...
            'censoring', Censored(Sample));
    
    %compute 'c' for bootstrap model on left out samples
    Cboot(i) = cIndex(Beta, X(:, Out).', Survival(Out), Censored(Out));
    
end
C = 0.368 * err + 0.632 * sum(Cboot) / B;
fprintf('\tTraining = %g, Bootstrap = %g, 0.632 c-index = %g\n',...
    err, sum(Cboot)/B, C);
