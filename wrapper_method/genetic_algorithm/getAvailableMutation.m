function y = getAvailableMutation(d)
% What will be done in this process:
%   1) Extract mutation data from training set
%   2) Drop all samples with number of NaN > 0


d.Features = d.Features(strcmp(d.SymbolTypes, 'Mutation'),strcmp(d.AvailableMutation, 'Yes'));
d.Symbols = d.Symbols(strcmp(d.SymbolTypes, 'Mutation'));
d.SymbolTypes = d.SymbolTypes(strcmp(d.SymbolTypes, 'Mutation'));

sample_keep_indices = (sum(isnan(d.Features), 1)==0);

d.Features=d.Features(:,sample_keep_indices);
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);
d.AvailableClinical = d.AvailableClinical(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);
d.AvailableMutation = d.AvailableMutation(sample_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(sample_keep_indices);

y = d;