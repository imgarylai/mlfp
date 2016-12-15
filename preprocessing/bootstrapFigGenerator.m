clear all
addpath ..
d=load('../data/BRCA.mat');
warning('off','all')
d.Symbols = cellstr(d.Symbols);
d.SymbolTypes = cellstr(d.SymbolTypes);
sample_keep_indices = ~isnan(d.Survival) & ~isnan(d.Censored);
d.Features = d.Features(:,sample_keep_indices); 
d.Survival = d.Survival(sample_keep_indices);
d.Censored = d.Censored(sample_keep_indices);
d.AvailableProtein = d.AvailableProtein(sample_keep_indices);
d.AvailableClinical = d.AvailableClinical(sample_keep_indices);
d.AvailableCNV = d.AvailableCNV(sample_keep_indices);
d.AvailableMutation = d.AvailableMutation(sample_keep_indices);
d.AvailablemRNA = d.AvailablemRNA(sample_keep_indices);

B = 100;

for j=137:236
    sample_keep_indices = ~isnan(d.Features(j,:));
    F = d.Features(:,sample_keep_indices); 
    S = d.Survival(sample_keep_indices);
    C = d.Censored(sample_keep_indices);
    Beta1 = coxphfit(F(j,:).', S, 'Censoring', C);
    err = cIndex(Beta1, F(j,:).', S, C);
    
    
    for i = 1:B
    %generate random sample index
   

    Sample = ceil(length(S) * rand(1,length(S)));
    while(sum(C(Sample)) == length(Sample))
        Sample = ceil(length(S) * rand(1,length(S)));
    end
    N = length(S);
    %generate out of bootstrap sample
    Out = setdiff([1:N], Sample);
    
    %build model
    Beta = coxphfit(F(j,Sample).', S(Sample),...
        'censoring', C(Sample));
    
    %compute 'c' for bootstrap model on left out samples
    Cboot(i) = cIndex(Beta, F(j,Out).', S(Out), C(Out));
    end
    feature_c(j) = 0.368 * err + 0.632 * sum(Cboot)/B;
end

 feature_keep_indices=feature_c>0.50;
 findex=feature_c<0.5;
 feature_c(findex)=1-feature_c(findex);
 plot(feature_c)
 
%  feature_keep_indices=feature_c>0.51;
%  Features=Features(feature_keep_indices,:);
%  SymbolTypes=SymbolTypes(feature_keep_indices);
%  Symbols=Symbols(feature_keep_indices);

  hold on
  for x=137:236
      if feature_c(x)>0.58
          txt=d.Symbols(x);
          text(x,feature_c(x),strrep(txt,'_','\_'));
      end
  end