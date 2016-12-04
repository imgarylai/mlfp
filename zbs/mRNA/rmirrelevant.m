function y = rmirrelevant(d)
warning('off','all')
B = 20;
for j=1:size(d.Features,1)
    for i = 1:B
    %generate random sample index
    F = d.Features(j,:);
    S = d.Survival;
    C = d.Censored;
    Sample = ceil(length(F) * rand(1,length(F)));
    while(sum(C(Sample)) == length(Sample))
        Sample = ceil(length(F) * rand(1,length(F)));
    end
    N = length(S);
    %generate out of bootstrap sample
    Out = setdiff([1:N], Sample);
    
    %build model
    Beta = coxphfit(F(Sample).', S(Sample),...
        'censoring', C(Sample));
    
    %compute 'c' for bootstrap model on left out samples
    Cboot(i,j) = cIndex(Beta, F(Out).', S(Out), C(Out));
    end
end
% 
 feature_c = sum(Cboot,1)/B;
 feature_keep_indices=feature_c>0.50;
 d.Features = d.Features(feature_keep_indices,:);
 d.Symbols = d.Symbols(feature_keep_indices);
 d.SymbolTypes=d.SymbolTypes(feature_keep_indices);

 y=d;
%  feature_keep_indices=feature_c>0.51;
%  Features=Features(feature_keep_indices,:);
%  SymbolTypes=SymbolTypes(feature_keep_indices);
%  Symbols=Symbols(feature_keep_indices);
%  hold on
%  for x=1:length(Symbols)
%      if feature_c(x)>0.56
%          txt=Symbols(x);
%          text(x,feature_c(x),strrep(txt,'_','\_'));
%      end
%  end













