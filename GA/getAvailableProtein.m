function y = getAvailableProtein(d)
% What will be done in this process:
%   1) Extract Protein data from training set 
%   2) Drop samples with number of NaN > 0
%   3) Drop feature with NaN > Threshold

d.Features = d.Features(strcmp(d.SymbolTypes, 'Protein'),strcmp(d.AvailableProtein, 'Yes'));
d.Symbols = d.Symbols(strcmp(d.SymbolTypes, 'Protein'));
d.SymbolTypes = d.SymbolTypes(strcmp(d.SymbolTypes, 'Protein'));


feature_keep_indices = (sum(isnan(d.Features),2)==0);
d.Features = d.Features(feature_keep_indices,:);
d.Symbols = d.Symbols(feature_keep_indices,:);
d.SymbolTypes = d.SymbolTypes(feature_keep_indices,:);



sample_keep_indices = (sum(isnan(d.Features),1)==0); 
d.Features=d.Features(:,sample_keep_indices);
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);
d.AvailableClinical = d.AvailableClinical(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);
d.AvailableMutation = d.AvailableMutation(sample_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(sample_keep_indices);

y = d;