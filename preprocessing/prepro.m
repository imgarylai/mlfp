function y = prepro(d)
% What will be done in this process: 
%   1) Drop all samples with unavailable Censored or Survivial
%   2) Convert Symbols and SymbolTypes to a readable format
%   3) Remove useless attributesv 'Samples'
%   4) Drop features that number of NaN greater than 500
%   5) 500 means we drop a feature if more than 500 samples do not have
%       that feature, 500 is a trade off value, choice 900 will have less
%       effect of dropping feature, choice 100 will have drop too many feature
%       due to the missing rate of AvailableXXX itself. 


d.Symbols = cellstr(d.Symbols);
d.SymbolTypes = cellstr(d.SymbolTypes);
feature_keep_indices = (sum(isnan(d.Features), 2) <500 & range(d.Features, 2)~= 0 );


sample_keep_indices = ~isnan(d.Survival) & ~isnan(d.Censored);
d.Features = d.Features(feature_keep_indices,sample_keep_indices); 
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);
d.AvailableClinical = d.AvailableClinical(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);
d.AvailableMutation = d.AvailableMutation(sample_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(sample_keep_indices);
d.SymbolTypes=d.SymbolTypes(feature_keep_indices);
d.Symbols=d.Symbols(feature_keep_indices);

y=d;



