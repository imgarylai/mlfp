function newFeatures = ml_PCA(mode, features_clinical, features_CNV, ... 
    features_mutation, features_mRNA_protein, Survival, Censored, selection_params) 
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
persistent features_clinical_PCA_mapping; 
persistent features_CNV_PCA_mapping; 
persistent features_mutation_PCA_mapping; 
persistent features_mRNA_protein_PCA_mapping; 



if (isempty(selection_params))
    N1 = 10; 
    N2 = 10; 
    N3 = 10;
    N4 = 10; 
else
    if (size(selection_params) ~= 4)
        error('we need four parameters for PCA');
    end 
    N1 = selection_params(1);
    N2 = selection_params(2);
    N3 = selection_params(3);
    N4 = selection_params(4);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PCA begins 
if (strcmp(mode, 'train'))
    [features_clinical_PCA, features_clinical_PCA_mapping] = ... 
        compute_mapping(features_clinical, 'PCA', N1); 
    [features_CNV_PCA, features_CNV_PCA_mapping] = ... 
        compute_mapping(features_CNV, 'PCA', N2); 
    [features_mutation_PCA, features_mutation_PCA_mapping] = ... 
        compute_mapping(features_mutation, 'PCA', N3); 
    [features_mRNA_protein_PCA, features_mRNA_protein_PCA_mapping] = ... 
        compute_mapping(features_mRNA_protein, 'PCA', N4); 
else
    % make sure it has been trained, 
    % otherwise you should not get mapped test features 
    if (isempty(features_clinical_PCA_mapping))
        error('you haven''t run training data through PCA yet, please do this first'); 
        
    end 
    features_clinical_PCA = features_clinical * features_clinical_PCA_mapping.M;
    features_clinical_PCA = features_clinical_PCA - mean(features_clinical_PCA); 
    features_CNV_PCA = features_CNV * features_CNV_PCA_mapping.M;
    features_CNV_PCA = features_CNV_PCA - mean(features_CNV_PCA); 
    features_mutation_PCA = features_mutation * features_mutation_PCA_mapping.M;
    features_mutation_PCA = features_mutation_PCA - mean(features_mutation_PCA); 
    features_mRNA_protein_PCA = features_mRNA_protein * features_mRNA_protein_PCA_mapping.M;
    features_mRNA_protein_PCA = features_mRNA_protein_PCA - mean(features_mRNA_protein_PCA); 
end 


newFeatures = [features_clinical_PCA features_CNV_PCA ... 
        features_mutation_PCA features_mRNA_protein_PCA]; 

