function y = rmirrelevant(d)
warning('off','all')

B = 30;

for j=1:size(d.Features,1)
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
        Sample = ceil(length(F) * rand(1,length(F)));
    end
    N = length(S);
    %generate out of bootstrap sample
    Out = setdiff([1:N], Sample);
    
    %build model
    Beta = coxphfit(F(Sample).', S(Sample),...
        'censoring', C(Sample));
    
    %compute 'c' for bootstrap model on left out samples
    Cboot(i) = cIndex(Beta, F(Out).', S(Out), C(Out));
    end
    feature_c(j) = 0.368 * err + 0.632 * sum(Cboot)/B;
end
% 

%  findex=feature_c<0.5;
%  feature_c(findex)=1-feature_c(findex);

 feature_keep_indices=feature_c>0.52;
 d.Features = d.Features(feature_keep_indices,:);
 d.Symbols = d.Symbols(feature_keep_indices);
 d.SymbolTypes=d.SymbolTypes(feature_keep_indices);

 y=d;













