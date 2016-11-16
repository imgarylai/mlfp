function y = preprocessing(d)
% What will be done in this process:
%   1) Drop all never changed features
%   2) Drop all samples with unavailable Censored or Survivial
%   3) Convert Symbols and SymbolTypes to a readable format
%   4) Remover Samples attributes (useless)

d.Symbols = cellstr(d.Symbols);
d.SymbolTypes = cellstr(d.SymbolTypes);
d = rmfield(d,'Samples');
% note that NaN will be automatic ignored when using range() function
feature_keep_indices = (range(d.Features, 2) ~= 0);
d.Features = d.Features(feature_keep_indices, :); 
d.Symbols = d.Symbols(feature_keep_indices);
d.SymbolTypes = d.SymbolTypes(feature_keep_indices);

sample_keep_indices = ~isnan(d.Survival) & ~isnan(d.Censored);
d.Features = d.Features(:,sample_keep_indices); 
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);
d.AvailableClinical = d.AvailableClinical(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);
d.AvailableMutation = d.AvailableMutation(sample_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(sample_keep_indices);
y=d;



