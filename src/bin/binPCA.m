

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% demonstration parameters 
addpath(genpath('../utils')); 
addpath(genpath('../feature_selections')); 


% func = @test_feature_selection; 
func = @ml_PCA; 

% evaluation_func = @eva; 
evaluation_func = @eva_cv; 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% method I I want to run evaluation myself and my func does not have params 
% for training and obtaining validation CI 
% disp('****************************************************************');
% disp('Method 1, no parameter is provided to feature selection function');
% CI_validation = evaluation_func(func, 1, []); 
% 
% % for testing 
% CI_test       = evaluation_func(func, 2, []);
% fprintf('validation CI %f, test CI %f\n', CI_validation, CI_test); 






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% method II I want to pass in a group of parameters and let the framework 
%% figures out what's the best parameter group 
% for training and obtaining validation CI 

N = 20; 
fprintf('\n')
disp('****************************************************************');
fprintf('PCA, %d groups of parameters are provided to feature selection function\n', N);
selection_params = zeros(N, 4); 
for i=1:N
    selection_params(i, :) = [i*5 i*5 i*5 i*5]; 
end 
    

[CI_validation, best_params] = evaluation_func(func, 1, selection_params);


% for testing
CI_test = evaluation_func(func, 2, best_params); 
fprintf('validation CI %f, test CI %f, best parameter ', CI_validation, CI_test); 
disp(best_params); 









