function [CI, best_params] = eva(feature_selection_func, eva_type, feature_select_params) 
% used for evaluate feature selection function 
% input: feature_selection_func: function pointer for feature selection 
%        eva_type: 1: train and validation 
%                  2: test 
%        feature_select_params: a matrix of size n*m, where n is the number
%        of parameter group, m is the number of parameters in each group 
% output: CI: CI 
%         best_params: the group of parameters provides best CI 


persistent need_to_load_data 
persistent features_clinical_train 
persistent features_CNV_train 
persistent features_mutation_train 
persistent features_mRNA_protein_train 
persistent survival_train
persistent censored_train

persistent features_clinical_validation 
persistent features_CNV_validation 
persistent features_mutation_validation 
persistent features_mRNA_protein_validation 
persistent survival_validation
persistent censored_validation

persistent features_clinical_test 
persistent features_CNV_test 
persistent features_mutation_test 
persistent features_mRNA_protein_test 
persistent survival_test
persistent censored_test 

persistent Beta 
% persistent persistent_best_CI
% persistent persistent_best_beta 
% persistent persistent_best_params
NO_OUTPUT = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if this is first time run the evaluation function, load data from disk 
if (isempty(need_to_load_data))
    load '../../data/clean/cleanedData_BRCA.mat'; 

    N = length(Survival); 


    % now seperate data into train, validation, test
    Train = 1 : ceil(N/2);
    Validation = ceil(N/2)+1 : ceil(3*N/4);
    Test = ceil(3*N/4)+1 : N; 

    features_clinical_train       =    features_clinical(Train, :); %#ok<*NODEF>
    features_CNV_train            =    features_CNV(Train, :);      %#ok<*NODEF>
    features_mutation_train       =    features_mutation(Train, :); %#ok<*NODEF>
    features_mRNA_protein_train   =    features_mRNA_protein(Train, :); %#ok<*NODEF>
    survival_train                =    Survival(Train); 
    censored_train                =    Censored(Train); 

    features_clinical_validation       =    features_clinical(Validation, :); %#ok<*NODEF>
    features_CNV_validation            =    features_CNV(Validation, :);      %#ok<*NODEF>
    features_mutation_validation       =    features_mutation(Validation, :); %#ok<*NODEF>
    features_mRNA_protein_validation   =    features_mRNA_protein(Validation, :); %#ok<*NODEF>
    survival_validation                =    Survival(Validation); 
    censored_validation                =    Censored(Validation); 

    features_clinical_test       =    features_clinical(Test, :); %#ok<*NODEF>
    features_CNV_test            =    features_CNV(Test, :);      %#ok<*NODEF>
    features_mutation_test       =    features_mutation(Test, :); %#ok<*NODEF>
    features_mRNA_protein_test   =    features_mRNA_protein(Test, :); %#ok<*NODEF>
    survival_test                =    Survival(Test); 
    censored_test                =    Censored(Test); 

    
    disp('load data') 
    need_to_load_data = 1;
end 



warning('off','all')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% begin feature selection test 
% feature_selection_func  =   @test_feature_selection;
if (eva_type == 1)
    % use train and validation to obatin the best beta 
    if ~isempty(feature_select_params) && size(feature_select_params, 1) ~= 0
        % if you choose to pass in parameters for feature selection func 
        CIs = zeros(size(feature_select_params, 1)); 
        for i=1:size(feature_select_params, 1)
            X_train  =  feature_selection_func('train', features_clinical_train, ...
                            features_CNV_train, features_mutation_train, ...
                            features_mRNA_protein_train, survival_train, ... 
                            censored_train, feature_select_params(i, :));   
            Beta = coxphfit(X_train, survival_train, 'Censoring', censored_train);


            X_validation = feature_selection_func('validation', features_clinical_validation, ...
                                features_CNV_validation, features_mutation_validation, ...
                                features_mRNA_protein_validation, survival_validation, ...
                                censored_validation, feature_select_params(i, :)); 

            CI_validation   =   cIndex(Beta, X_validation, survival_validation, censored_validation);
            if ~NO_OUTPUT
                fprintf('group %d of parameters, standard validation CI: %f\n', i, CI_validation); 
            end 
            CIs(i) = CI_validation;             
        end
        
        % now select the best CI and corresponding parameter group 
        [max_CI, max_CI_index] = max(CIs); 
        max_CI = max_CI(1); 
        max_CI_index = max_CI_index(1); 
        CI = max_CI; 
        best_params = feature_select_params(max_CI_index, :); 
        fprintf('group %d has the best performance, standard validation CI:%f\n', max_CI_index, max_CI); 
        
    % no feature selection parameter is provided 
    else
        X_train  =  feature_selection_func('train', features_clinical_train, ...
                        features_CNV_train, features_mutation_train, ...
                        features_mRNA_protein_train, survival_train, censored_train, []);   
        Beta = coxphfit(X_train, survival_train, 'Censoring', censored_train);


        X_validation = feature_selection_func('validation', features_clinical_validation, ...
                            features_CNV_validation, features_mutation_validation, ...
                            features_mRNA_protein_validation, survival_validation, ...
                            censored_validation, []); 

        CI_validation   =   cIndex(Beta, X_validation, survival_validation, censored_validation);
        CI = CI_validation; 
        best_params = []; 
        if ~NO_OUTPUT
            fprintf('standard validation CI: %f\n', CI_validation); 
        end 
    end
elseif (eva_type == 2)
    % test 
    if (~isempty(feature_select_params) && size(feature_select_params, 1) ~= 0)
        if (size(feature_select_params, 1) ~= 1)
            error('for test only one group parameter is allowed');
        end
        X_test  =  feature_selection_func('test', features_clinical_test, ...
                features_CNV_test, features_mutation_test, ...
                features_mRNA_protein_test, survival_test, censored_test, feature_select_params);         
    else
        X_test  =  feature_selection_func('test', features_clinical_test, ...
                features_CNV_test, features_mutation_test, ...
                features_mRNA_protein_test, survival_test, censored_test, []); 
    end 
    
    CI_test =  cIndex(Beta, X_test, survival_test, censored_test); 
    CI = CI_test; 
    best_params = []; 
    if ~NO_OUTPUT
        fprintf('standard test CI: %f\n', CI); 
    end
else
    disp('what do you give as type');
end

end


