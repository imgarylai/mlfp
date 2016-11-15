

% this file is a demonstration of how to use the framework to evaluate your
% feature selection function 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% guide on how to write your own feature selection func 
% please refer to test_feature_selection.m and ml_PCA.m 





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demonstration parameters 
addpath(genpath('utils')); 
addpath(genpath('feature_selections')); 


func = @test_feature_selection; 
% func = @ml_PCA; 

evaluation_func = @eva; 
% evaluation_func = @eva_cv; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% method I I want to run evaluation myself and my func does not have params 
% for training and obtaining validation CI 
disp('****************************************************************');
disp('Method 1, no parameter is provided to feature selection function');
CI_validation = evaluation_func(func, 1, []); 

% for testing 
CI_test       = evaluation_func(func, 2, []);
fprintf('validation CI %f, test CI %f\n', CI_validation, CI_test); 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% method II I want to pass in a group of parameters and let the framework 
%% figures out what's the best parameter group 
% for training and obtaining validation CI 

fprintf('\n')
disp('****************************************************************');
disp('Method 2, three groups of parameters are provided to feature selection function');
if (isequal(func, @test_feature_selection))
    [CI_validation, best_params] = evaluation_func(func, 1, [10; 20; 30]);
elseif (isequal(func, @ml_PCA)) 
    selection_params = [10 10 10 10; 20 20 20 20; 30 30 30 30]; 
    [CI_validation, best_params] = evaluation_func(func, 1, selection_params);
end 


% for testing
% we need to run the framework one more time to get the best beta 
% CI_validation = evaluation_func(func, 1, best_params); 
CI_test = evaluation_func(func, 2, best_params); 
fprintf('validation CI %f, test CI %f, best parameter ', CI_validation, CI_test); 
disp(best_params); 









