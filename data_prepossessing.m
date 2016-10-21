

%load data 
clc; clear; 
load BRCA.Data.mat; 
maxNumCompThreads(8); 
warning('off','all')
Symbols = cellstr(Symbols);
SymbolTypes = cellstr(SymbolTypes); 


fprintf('number of Nan in Censored %d\n', sum(isnan(Censored)));
fprintf('number of Nan in Survival %d\n\n', sum(isnan(Survival)));

% within a feature, if number of missing(Nan) larger than 
% total records * DROP_FEATURE_THRESHOLD, then drop current feature 
fprintf('before dropping features with too many Nan, number of features %d\n', size(Features, 1))
DROP_FEATURE_THRESHOLD = 1/10; 
feature_keep_indices = (sum(isnan(Features), 2)) < ...
    length(Censored) * DROP_FEATURE_THRESHOLD; 
features_keep = Features(feature_keep_indices, :); 
Symbols = Symbols(feature_keep_indices);
SymbolTypes = SymbolTypes(feature_keep_indices);
fprintf('after dropping features with too many Nan, number of features %d\n\n', size(features_keep, 1))

% now drop all features with range == 0 
fprintf('before dropping features with range 0, number of features %d\n', size(features_keep, 1))
feature_keep_indices = (range(features_keep, 2) ~= 0);
features_keep = features_keep(feature_keep_indices, :); 
Symbols = Symbols(feature_keep_indices);
SymbolTypes = SymbolTypes(feature_keep_indices);
fprintf('after dropping features with range 0, number of features %d\n\n', size(features_keep, 1))


% imputation using euclidean distance with k neighbours 
fprintf('before imputation, number of Nan %d\n', sum(sum(isnan(features_keep))));
IMPUTATION_K = 1; 
% features_imputed = knnimpute(features_keep', IMPUTATION_K)'; 
fprintf('after imputation, number of Nan %d\n\n', sum(sum(isnan(features_keep))));

% save('new_data_BRCA.mat')


% begin feature selection 

M = 200;
X = features_keep(1:M, :); 


Keep = ~isnan(Survival) & ~isnan(Censored) & (sum(isnan(X), 1) == 0);
X = X(:, Keep);
Survival = Survival(Keep);
Censored = Censored(Keep);



save('newData.mat')
