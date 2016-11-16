function [CI, best_params] = eva_cv(feature_selection_func, eva_type, feature_select_params)
% 3 fold cv used for evaluate feature selection function
% input: feature_selection_func: function pointer for feature selection
%        eva_type: 1: train and validation
%                  2: test
%        feature_select_params: a matrix of size n*m, where n is the number
%        of parameter group, m is the number of parameters in each group
% output: CI: CI
%         best_params: the group of parameters provides best CI


persistent need_to_load_data
persistent features_clinical_train_val
persistent features_CNV_train_val
persistent features_mutation_train_val
persistent features_mRNA_protein_train_val
persistent survival_train_val
persistent censored_train_val


persistent features_clinical_test
persistent features_CNV_test
persistent features_mutation_test
persistent features_mRNA_protein_test
persistent survival_test
persistent censored_test


persistent Folds;


persistent Beta
NO_OUTPUT = false;
K = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% if this is first time run the evaluation function, load data from disk
if (isempty(need_to_load_data))
    load 'data/clean/cleanedData_BRCA.mat'; 


    N = length(Survival);
    N_Train_Validation = ceil(N / 4 * 3);
    Train_Validation = 1 : ceil(N / 4 * 3);
    Test = ceil(3*N/4)+1 : N;


    Folds = ceil(Train_Validation / (N_Train_Validation/K));
    Folds(Folds>K) = 0;



    features_clinical_train_val       =    features_clinical(Train_Validation, :); %#ok<*NODEF>
    features_CNV_train_val            =    features_CNV(Train_Validation, :);      %#ok<*NODEF>
    features_mutation_train_val       =    features_mutation(Train_Validation, :); %#ok<*NODEF>
    features_mRNA_protein_train_val   =    features_mRNA_protein(Train_Validation, :); %#ok<*NODEF>
    survival_train_val                =    Survival(Train_Validation);
    censored_train_val                =    Censored(Train_Validation);


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
if (eva_type == 1)
    % use train and validation to obatin the best beta
    if ~isempty(feature_select_params) && size(feature_select_params, 1) ~= 0
        % if you choose to pass in parameters for feature selection func
        CIs = zeros(size(feature_select_params, 1));
        CI_stds = zeros(size(feature_select_params, 1));
        for i=1:size(feature_select_params, 1)
            %cycle through folds
            innerCIs = nan(1, K);
            for j = 1:K
                X_train  =  feature_selection_func('train', ...
                                features_clinical_train_val(Folds ~= j, :), ...
                                features_CNV_train_val(Folds ~= j, :), ...
                                features_mutation_train_val(Folds ~= j, :), ...
                                features_mRNA_protein_train_val(Folds ~= j, :), ...
                                survival_train_val(Folds ~= j, :), ...
                                censored_train_val(Folds ~= j, :), ...
                                feature_select_params(i, :));
                Beta = coxphfit(X_train, survival_train_val(Folds ~= j, :), ...
                    'Censoring', censored_train_val(Folds ~= j, :));


                X_validation = feature_selection_func('validation', ...
                                features_clinical_train_val(Folds == j, :), ...
                                features_CNV_train_val(Folds == j, :), ...
                                features_mutation_train_val(Folds == j, :), ...
                                features_mRNA_protein_train_val(Folds == j, :), ...
                                survival_train_val(Folds == j, :), ...
                                censored_train_val(Folds == j, :), ...
                                feature_select_params(i, :));

                CI_validation  =  cIndex(Beta, X_validation, ...
                                    survival_train_val(Folds == j, :), ...
                                    censored_train_val(Folds == j, :));
                innerCIs(j) = CI_validation;
            end

            CI          =       mean(innerCIs);
            CI_std      =       std(innerCIs);
            CIs(i)      =       CI;
            CI_stds(i)  =       CI_std;

            if ~NO_OUTPUT
                fprintf('>>>group %d of parameters, %d fold cv CI: %f, std %f\n', ...
                            i, K, CI, CI_std);
            end
        end

        % now select the best CI and corresponding parameter group
        [max_CI, max_CI_index] = max(CIs);
        max_CI = max_CI(1);
        max_CI_index = max_CI_index(1);
        CI = max_CI;
        CI_std = CI_stds(max_CI_index);
        best_params = feature_select_params(max_CI_index, :);
        fprintf('>>>group %d has the best performance, %d fold cross validation CI:%f, std %f\n', ...
                    max_CI_index, K, CI, CI_std);

    % no feature selection parameter is provided
    else
        %cycle through folds
        innerCIs = nan(1, K);
        for i = 1:K
            X_train  =  feature_selection_func('train', ...
                            features_clinical_train_val(Folds ~= i, :), ...
                            features_CNV_train_val(Folds ~= i, :), ...
                            features_mutation_train_val(Folds ~= i, :), ...
                            features_mRNA_protein_train_val(Folds ~= i, :), ...
                            survival_train_val(Folds ~= i, :), ...
                            censored_train_val(Folds ~= i, :), []);
            Beta = coxphfit(X_train, survival_train_val(Folds ~= i, :), ...
                'Censoring', censored_train_val(Folds ~= i, :));


            X_validation = feature_selection_func('validation', ...
                            features_clinical_train_val(Folds == i, :), ...
                            features_CNV_train_val(Folds == i, :), ...
                            features_mutation_train_val(Folds == i, :), ...
                            features_mRNA_protein_train_val(Folds == i, :), ...
                            survival_train_val(Folds == i, :), ...
                            censored_train_val(Folds == i, :), []);

            CI_validation  =  cIndex(Beta, X_validation, ...
                            survival_train_val(Folds == i, :), ...
                            censored_train_val(Folds == i, :));
            innerCIs(i) = CI_validation;
        end
        CI = mean(innerCIs);
        CI_std = std(innerCIs);

        best_params = [];
        if ~NO_OUTPUT
            fprintf('>>>%d fold cv CI: %f, std %f\n', K, CI, CI_std);
        end
    end
elseif (eva_type == 2)
    % test
    if (isempty(feature_select_params) || size(feature_select_params, 1) == 0)
        feature_select_params = [];
    else
        if (size(feature_select_params, 1) ~= 1)
            error('for test only one group parameter is allowed');
        end
    end

    % we need to retrain the model using all training data
    X_train  =  feature_selection_func('train', ...
                    features_clinical_train_val, ...
                    features_CNV_train_val, ...
                    features_mutation_train_val, ...
                    features_mRNA_protein_train_val, ...
                    survival_train_val, ...
                    censored_train_val, feature_select_params);
    Beta = coxphfit(X_train, survival_train_val, ...
        'Censoring', censored_train_val);

    % test
    X_test  =  feature_selection_func('test', features_clinical_test, ...
            features_CNV_test, features_mutation_test, ...
            features_mRNA_protein_test, survival_test, censored_test, feature_select_params);

    CI_test =  cIndex(Beta, X_test, survival_test, censored_test);
    CI = CI_test;
    best_params = feature_select_params;
    if ~NO_OUTPUT
        fprintf('>>>%d fold test CI: %f\n', K, CI);
    end
else
    disp('what do you give as type');
end

end
