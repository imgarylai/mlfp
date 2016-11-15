function newFeatures = test_feature_selection(mode, features_clinical, features_CNV, ... 
    features_mutation, features_mRNA_protein, Survival, Censored, selection_params) 
% this is an example function demonstrates how to select/extract 
% features and how to use the evaluation framework, it just select first N
% features_clinical parameters, default 20  
% REMEMBER, 
% you should only use training set to select/extract features, 
% don't touch validation and test set at this time, 
% in other words, your function should behave DIFFERENTLY 
% according different modes (train/validation, test) 
%
% input: mode: train/validation/test 
% features: each row is a sample, each column is a feature;
% Survival, Censored: 
% the reason Survival and Censored are passed in as parameters is that for
% some supervised feature extraction method, they are needed, but DON'T 
% use these two when you are dealing with validation and test data set  


% this test functions select first M features from all features without
% Nan, and return them as new features 

if (isempty(selection_params))
    M = 20;
else
    M = selection_params(1);
end
newFeatures = features_clinical(:, 1:M); 


%end 