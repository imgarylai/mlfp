%
% Input matrix: X, which is samples by features
%       
function [newFeatures, mseError]= ml_autoencoder(X, N)
    [originalFeatureSize, originalSampleSize] = size(X);
    disp('Please check the following information.');
    fprintf('Feature Number: %d, Sample Number: %d. \n', ...
        originalFeatureSize, originalSampleSize);
    fprintf('Reduce your matrix to %d dimension. \n', N);
    disp('Start training...');
    autoenc = trainAutoencoder(X, 'hiddenSize', N, 'MaxEpochs', 1000, ...
        'ShowProgressWindow', false);
    disp('Training done!');
    newFeatures = encode(autoenc, X);
    [newFeatureSize, newSampleSize] = size(newFeatures);
    fprintf('Feature Number: %d, Sample Number: %d. \n', newFeatureSize, newSampleSize);
    XReconstructed = predict(autoenc, X);
    mseError = mse(X - XReconstructed);
    fprintf('MSE error of Autoencoder: %f. \n', mseError);
%     v.features = newFeatures;
%     v.autoencoder = autoenc;
%     fname = sprintf('data/save/autoencoder-%d-%d.mat', originalFeatureSize, newFeatureSize);
%     save(fname, '-struct', 'v'); 
end