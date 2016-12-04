function [newFeatures,W1,b1] = sparse_autoencoder(input,hiddenSize)
    
    [originalFeatureSize, originalSampleSize] = size(input);
    sparsityParam = 0.01;
    lambda = 0.0001;
    beta = 3;

    disp('Please check the following information.');
    fprintf('Feature Number: %d, Sample Number: %d,Sparsity Parameter: %d,lambda: %d, beta: %d. \n', ...
        originalFeatureSize, originalSampleSize,sparsityParam,lambda,beta);
    fprintf('Reduce your matrix to %d dimension. \n', hiddenSize);

    disp('Start training...');
    [newFeatures,cost,W1,b1] = train(originalFeatureSize,hiddenSize,sparsityParam,lambda,beta,input);    
    disp('Training done!');
    [newFeatureSize, newSampleSize] = size(newFeatures);
    fprintf('cost of Autoencoder: %f. \n', cost);
end