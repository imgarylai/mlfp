function y = c_index_fitness(pop,d)
warning('off','all')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5-fold cross validation of a basic model (no model selection) %%%%%%%%%%
N = length(d.F);
Basic =d.F;
FeatIndex = find(pop==1); %Feature Index
Basic = Basic([FeatIndex],:);
%define number of folds and assign samples
K = 5;
Folds = ceil([1:N] / (N/K));
%cycle through folds, training and testing model
C = nan(1,K);
for i = 1:K
    Beta = coxphfit(Basic(:, Folds ~= i).', d.S(Folds ~= i).',...
        'Censoring', d.C(Folds ~= i).');
    C(i) = cIndex(Beta, Basic(:, Folds == i).', d.S(Folds == i),...
        d.C(Folds == i)); 
end

y=-mean(C);

