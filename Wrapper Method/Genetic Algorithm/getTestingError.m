function y = getTestingError(train,test,selected_feature)
% This function is used to report final testing error
% Here all seperated model will be combined to form a complete model
%
%
%
% remaining unfinished!!!!!
%
%
%

BasicTrain = train.Features([selected_feature],:);

BasicTest = test.Features([selected_feature],:);

Beta_Final = coxphfit(BasicTrain(:, :).', train.Survival(:).',...
        'Censoring', train.Censored(:).');
c_Final = cIndex(Beta_Final, BasicTest(:, :).', test.Survival(:),...
    test.Censored(:));

y=c_Final;
