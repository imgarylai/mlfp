

clc; clear; 
load './newData.mat'

N = length(Survival); 


% filter our features without 0 
Y = X~=0; 
feature_keep_indices = ~all(Y, 2); 

fprintf('before filtering features without 0, there are %d features\n', size(X, 1));
X = X(feature_keep_indices, :); 
fprintf('after  filtering features without 0, there are %d features\n', size(X, 1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Standard train/test/validation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Model assessment and selection with train/test/validation approach.\n');

%define train, test, and validation sets
Train = [1 : ceil(N/2)];
Validation = [ceil(N/2)+1 : ceil(3*N/4)];
Test = [ceil(3*N/4)+1 : N];


%build each model on training
Beta1 = coxphfit(X(:, Train).', Survival(Train).',...
                 'Censoring', Censored(Train).');

%compare performance of models on validation
C1 = cIndex(Beta1, X(:, Validation).', Survival(Validation),...
    Censored(Validation));

%report performance of best model on validation
    C = cIndex(Beta1, X(:, Test).', Survival(Test), Censored(Test));
    fprintf('\tModel 1 - without NOTCH1: c-index = %g\n', C);

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