function [newFeatures, mapping] = newPCA(features, selection_params)
% this is an example function that using PCA to show how to select/extract
% features and how to use the evaluation framework,
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
% selection_params: a vector contains all the parameters for selection



addpath(genpath('../extensions/drtoolbox'));



desired_dimentionality = selection_params; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA begins
% PCA uses the features passed in during training to find out the mapping from
% old features to new features, then the mapping is stored in a persistent
% variable features_PCA_mapping, then when the mode is test or validation,
% it just uses the variable features_PCA_mapping to transform features 


[newFeatures, features_PCA_mapping] = ...
            compute_mapping(features, 'PCA', desired_dimentionality);


mapping = features_PCA_mapping; 


    
% for test 
newFeatures = features * mapping.M; 
newFeatures = newFeatures - repmat(mean(newFeatures, 1), size(newFeatures, 1), 1);    

end
