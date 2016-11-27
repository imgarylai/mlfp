load 'data/clean/clean.mat';
T = 10;
X = clinical_train.Features;
[featureSize, sampleSize] = size(X);

for t = 1:T
    for i=1:featureSize
        [newFeatures, autoenc, mseError]= ml_autoencoder(X, i);
        if(t == 1)
            [newFeatureSize, newSampleSize] = size(newFeatures);
            v.features = newFeatures;
            v.autoencoder = autoenc;
            fname = sprintf('auto_encoder/generated_data/autoencoder-%s-%d-%d.mat', char(clinical_train.SymbolTypes(1)) ,featureSize, newFeatureSize);
            save(fname, '-struct', 'v'); 
        end
        
        mse(t, i) = mseError;
    end
end
