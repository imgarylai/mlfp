function y = getAvailableCNV(d)
% What will be done in this process:
%   1) Extract CNV data from training set 
%   2) Drop all samples with number of NaN > 0
d = rmfield(d,'AvailableProtein');
d = rmfield(d,'AvailableClinical');
d = rmfield(d,'AvailablemRNA');
d = rmfield(d,'AvailableMutation');
d = rmfield(d,'Samples');

d.Features = d.Features(strmatch({'CNV'},d.SymbolTypes),strcmp(d.AvailableCNV, 'Yes'));
d.Symbols = d.Symbols(strmatch({'CNV'},d.SymbolTypes));
d.SymbolTypes = d.SymbolTypes(strmatch({'CNV'},d.SymbolTypes));

% 2 is a tradeoff value for increasing many sample without lossing many feature
sample_keep_indices = (sum(isnan(d.Features), 1)<10);
d.Features=d.Features(:,sample_keep_indices);
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);

feature_keep_indices= (sum(isnan(d.Features), 2)==0 & range(d.Features,2)~=0);
d.Features=d.Features(feature_keep_indices,:);
d.SymbolTypes=d.SymbolTypes(feature_keep_indices);
d.Symbols=d.Symbols(feature_keep_indices);
y = d;



