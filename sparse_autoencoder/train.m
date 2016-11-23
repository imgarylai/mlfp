clear;
visibleSize     =  4;   
hiddenSize      =  2 ;  
sparsityParam   = 0.8 ; % desired average activation of the hidden units.                     
%lambda      = 0.0001;    % weight decay parameter       
%beta        = 9;         % weight of sparsity penalty term       
input = [0,0,0,1;0,0,1,0;0,1,0,0;1,0,0,0];
record = [];
opttheta = [];
optbeta = 0;
mincost = 9999999;

for sparsityParam = 0:0.1:1
    for lambda = 0:0.1:1
        for beta = 1:0.1:10
            %initialize the parameters
            theta = initializeParameters(hiddenSize, visibleSize);

            %  Use minFunc to minimize the function
            addpath minFunc/
            options.Method = 'lbfgs';
            options.maxIter = 400;	   
            options.display = 'on';
            [opttheta, cost] = minFunc( @(p) sparseAutoencoderCost(p,visibleSize, hiddenSize,lambda, sparsityParam,beta, input),theta, options);
            if cost < mincost
                mincost = cost;
                opttheta_best = opttheta;
                optbeta = beta;
                optlambda = lambda;
                optsparsityParam = sparsityParam;

            end
        end
    end
end

%[a b] = min(record);
%opttheta_best = record(b(2),1);
%optbeta = record(b(2),3);

W1 = reshape(opttheta(1:hiddenSize * visibleSize), hiddenSize, visibleSize);
b1 = opttheta(2*hiddenSize*visibleSize+1:2*hiddenSize*visibleSize+hiddenSize);
z2 = W1 * input + repmat(b1,[1, size(input,2)]);

answer = sigmoid(z2);