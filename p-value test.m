% Without normalization
Keep1 = ~isnan(Survival);
Keep2 = ~isnan(Censored);
Keep=Keep1 & Keep2;
F= Features(:,Keep)';
S= Survival(:,Keep)';
C= Censored(:,Keep)';
for i=1:17584
    [b,logl,H,stats] = coxphfit(F(:,i),S(:,1),'Censoring',C(:,1));
    P(1,i)=stats.p;
end
sum(P(1,:)<0.05)
% 1904
