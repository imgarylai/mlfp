function y = getAvailableProtein(d)
% What will be done in this process:
%   1) Extract clinical data from training set 
%   2) Drop all samples with number of NaN > 0
d = rmfield(d,'AvailableCNV');
d = rmfield(d,'AvailableMutation');
d = rmfield(d,'AvailablemRNA');
d = rmfield(d,'AvailableClinical');

d.Features = d.Features(strcmp(d.SymbolTypes, 'Protein'),strcmp(d.AvailableProtein, 'Yes'));
d.Symbols = d.Symbols(strcmp(d.SymbolTypes, 'Protein'));
d.SymbolTypes = d.SymbolTypes(strcmp(d.SymbolTypes, 'Protein'));

% 2 is a tradeoff value for increasing many sample without lossing many feature
sample_keep_indices = (sum(isnan(d.Features), 1)<10);
d.Features=d.Features(:,sample_keep_indices);
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);

feature_keep_indices= (sum(isnan(d.Features), 2)==0 & range(d.Features,2)~=0);
d.Features=d.Features(feature_keep_indices,:);
d.SymbolTypes=d.SymbolTypes(feature_keep_indices);
d.Symbols=d.Symbols(feature_keep_indices);

y = d;



