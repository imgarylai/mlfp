load('cleanedData_BRCA.mat');
%change the dataset here
dataset = features_clinical;
%use lasso to fit
[B,FitInfo] = lasso(dataset, Survival,'CV',10);
beltahat = B(:,FitInfo.IndexMinMSE);
outputfeature = [];
[H,W] = size(dataset);
%pick up the feature
for i = 1:W
    if beltahat(i,:) ~= 0
        outputfeature = [outputfeature,dataset(:,i)];
    end
end