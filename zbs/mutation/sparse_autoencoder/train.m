function [answer,cost,W1,b1] = train(visibleSize,hiddenSize,sparsityParam,lambda,beta,input)    
    %visibleSize =  4;   
    %hiddenSize =  2 ;  
    %sparsityParam = 0.01 ; % desired average activation of the hidden units.                     
    %lambda = 0.0001;    % weight decay parameter       
    %beta = 3;         % weight of sparsity penalty term       
    %input = [0,0,0,1;0,0,1,0;0,1,0,0;1,0,0,0];
    
    % initialize the parameters
    theta = initializeParameters(hiddenSize, visibleSize);

    % sue minFunc to minimize the function
    addpath minFunc/
    options.Method = 'lbfgs';
    options.maxIter = 400;	   
    options.display = 'on';
    [opttheta, cost] = minFunc( @(p) sparseAutoencoderCost(p,visibleSize, hiddenSize,lambda, sparsityParam,beta, input),theta, options);

    % get the hidden layer
    W1 = reshape(opttheta(1:hiddenSize * visibleSize), hiddenSize, visibleSize);
    b1 = opttheta(2*hiddenSize*visibleSize+1:2*hiddenSize*visibleSize+hiddenSize);
    z2 = W1 * input + repmat(b1,[1, size(input,2)]);
    answer = sigmoid(z2);

end
