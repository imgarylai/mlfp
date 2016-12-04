
visibleSize     =  4;  % number of input units 
hiddenSize      =  2 ;  % number of hidden units 
sparsityParam   = 0.01 ; % desired average activation of the hidden units.
                         % (This was denoted by the Greek alphabet rho, which looks like a lower-case "p",
                         %  in the lecture notes). 
lambda      = 0.0001;    % weight decay parameter       
beta        = 3;         % weight of sparsity penalty term       
%numpatches  = 10000;

%%======================================================================

%patches = sampleIMAGES(numpatches);
%[h, arr] = display_network(patches(:,randi(size(patches,2),208,1)),8);
input = [0,0,0,1;0,0,1,0;0,1,0,0;1,0,0,0];

%  Obtain random parameters theta
theta = initializeParameters(hiddenSize, visibleSize);

%%======================================================================

[cost, grad] = sparseAutoencoderCost(theta, visibleSize, hiddenSize, lambda, ...
                                     sparsityParam, beta, input);
%%======================================================================

checkNumericalGradient(); 
numgrad = computeNumericalGradient( @(x) sparseAutoencoderCost(x, visibleSize, ...
                                                  hiddenSize, lambda, ...
                                                  sparsityParam, beta, ...
                                                  input), theta);
diff = norm(numgrad-grad)/norm(numgrad+grad);
disp(diff); % Should be small. In our implementation, these values are
            % usually less than 1e-9.

            % When you got this working, Congratulations!!! 

