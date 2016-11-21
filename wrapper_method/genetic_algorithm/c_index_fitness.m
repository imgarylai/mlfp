function y = c_index_fitness(pop,d)
warning('off','all')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5-fold cross validation of a basic model (no model selection) %%%%%%%%%%
N = length(d.Features);
Basic =d.Features;
FeatIndex = find(pop==1); %Feature Index
Basic = Basic([FeatIndex],:);
%define number of folds and assign samples
K = 5;
Folds = ceil([1:N] / (N/K));
%cycle through folds, training and testing model
C = nan(1,K);
for i = 1:K
    Beta = coxphfit(Basic(:, Folds ~= i).', d.Survival(Folds ~= i).',...
        'Censoring', d.Censored(Folds ~= i).');
    C(i) = cIndex(Beta, Basic(:, Folds == i).', d.Survival(Folds == i),...
        d.Censored(Folds == i));
end
y=-mean(C);
