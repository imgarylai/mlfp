function newFeatures = ml_PCA(mode, features, feature_type, selection_params)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% parameters for PCA 
% this is used to store the the weight or mapping parameter obtained from 
% training,
% a persistent variable will keep its value, after you finish this
% function, in other words, declare all the parameters/weights that
% obtained during training as persistent, then you are able to use it
% during validation/test 

persistent features_clinical_PCA_mapping;
persistent features_CNV_PCA_mapping;
persistent features_mutation_PCA_mapping;
persistent features_mRNA_protein_PCA_mapping;



% this is the passed in parameter for PCA, 
% in PCA this parameter is just a scalar that specifies the desired
% dimentionality of newX 

desired_dimentionality = selection_params; 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA begins
% PCA uses the features passed in during training to find out the mapping from
% old features to new features, then the mapping is stored in a persistent
% variable features_PCA_mapping, then when the mode is test or validation,
% it just uses the variable features_PCA_mapping to transform features 

if (strcmp(mode, 'train'))
    if (strcmp(feature_type, 'clinical'))
        [newFeatures, features_clinical_PCA_mapping] = ...
            compute_mapping(features, 'PCA', desired_dimentionality);
    elseif (strcmp(feature_type, 'CNV'))
        [newFeatures, features_CNV_PCA_mapping] = ...
            compute_mapping(features, 'PCA', desired_dimentionality);
    elseif (strcmp(feature_type, 'mutation'))
        [newFeatures, features_mutation_PCA_mapping] = ...
            compute_mapping(features, 'PCA', desired_dimentionality);
    elseif (strcmp(feature_type, 'mRNA_protein'))
        [newFeatures, features_mRNA_protein_PCA_mapping] = ...
            compute_mapping(features, 'PCA', desired_dimentionality);            
    end
else
    % make sure it has been trained, meaning features_mapping is not
    % uninitialized 
    % otherwise you should not get mapped test features

    %%%%%% select corresponding mapping according to feature type 
    if (strcmp(feature_type, 'clinical'))
        current_mapping = features_clinical_PCA_mapping; 
    elseif (strcmp(feature_type, 'CNV'))
        current_mapping = features_CNV_PCA_mapping; 
    elseif (strcmp(feature_type, 'mutation'))
        current_mapping = features_mutation_PCA_mapping; 
    elseif (strcmp(feature_type, 'mRNA_protein'))
        current_mapping = features_mRNA_protein_PCA_mapping; 
    end

    if (isempty(current_mapping))
        error('you haven''t run training data through PCA yet, please do this first');
    end
    
    newFeatures = features * current_mapping.M; 
    newFeatures = newFeatures - repmat(mean(newFeatures, 1), size(newFeatures, 1), 1);    
end
