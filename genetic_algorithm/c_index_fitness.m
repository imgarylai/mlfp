function y = c_index_fitness(pop,Train)
warning('off','all')
addpath ..
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 5-fold cross validation of a basic model (no model selection) %%%%%%%%%%
% N = length(d.F);
% Basic =d.F;
 FeatIndex = find(pop==1); %Feature Index
 Train.F = Train.F([FeatIndex],:);
  
 
 
% %define number of folds and assign samples
% K = 5;
% Folds = ceil([1:N] / (N/K));
% %cycle through folds, training and testing model
% C = nan(1,K);
% for i = 1:K


Beta = coxphfit(Train.F.', Train.S.',...
        'Censoring', Train.C.');
C = cIndex(Beta,  Train.F.', Train.S,...
         Train.C); 
% end

y=-C;
