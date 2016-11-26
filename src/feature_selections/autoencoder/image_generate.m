load 'data/clean/clean.mat';
T = 10;
X = clinical_train.Features;
[featureSize, sampleSize] = size(X);

for t = 1:T
    mse = zeros(featureSize, 1);
    for i=1:featureSize
        [newFeatures, mseError]= ml_autoencoder(X, i);
        mse(t, i) = mseError;
    end
end