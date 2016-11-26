function [newFeatures,cost] = ml_sparse_autoencoder(input,hiddenSize,sparsityParam,lambda,beta)
    [originalFeatureSize, originalSampleSize] = size(input);
    %sparsityParam = 0.01;
    %lambda = 0.0001;
    %beta = 3;
    
   
    disp('Please check the following information.');
    fprintf('Feature Number: %d, Sample Number: %d,Sparsity Parameter: %d,lambda: %d, beta: %d. \n', ...
        originalFeatureSize, originalSampleSize,sparsityParam,lambda,beta);
    fprintf('Reduce your matrix to %d dimension. \n', hiddenSize);
    disp('Start training...');
    [newFeatures,cost,W1,b1] = train(originalFeatureSize,hiddenSize,sparsityParam,lambda,beta,input);    
    disp('Training done!');
    [newFeatureSize, newSampleSize] = size(newFeatures);
    fprintf('Feature Number: %d, Sample Number: %d. \n', newFeatureSize, newSampleSize);
    fprintf('cost of Autoencoder: %f. \n', cost);
    v.features = newFeatures;
    v.weight = W1;
    v.b = b1;
    fname = sprintf('data/save/autoencoder-%d-%d.mat', originalFeatureSize, newFeatureSize);
    save(fname, '-struct', 'v');

end

