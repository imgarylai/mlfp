function [cost,grad] = sparseAutoencoderCost(theta, visibleSize, hiddenSize,lambda, sparsityParam, beta, data)

W1 = reshape(theta(1:hiddenSize*visibleSize), hiddenSize, visibleSize);
W2 = reshape(theta(hiddenSize*visibleSize+1:2*hiddenSize*visibleSize), visibleSize, hiddenSize);
b1 = theta(2*hiddenSize*visibleSize+1:2*hiddenSize*visibleSize+hiddenSize);
b2 = theta(2*hiddenSize*visibleSize+hiddenSize+1:end);

cost = 0;
W1grad = zeros(size(W1)); 
W2grad = zeros(size(W2));
b1grad = zeros(size(b1)); 
b2grad = zeros(size(b2));

m = size(data,2);

a1 = data;
z2 = W1 * data + repmat(b1,[1, size(data,2)]);
a2 =  sigmoid(z2);   
z3 = W2 * a2 + repmat(b2,[1, size(a2,2)]);
a3 =  sigmoid(z3);

% sparsity constratints 
rho_hat = mean(a2,2);                               
%hidden-layer-size x 1 
sparsity_term = - sparsityParam./rho_hat + (1-sparsityParam)./(1-rho_hat);

% error matrix
delta_3 = - (data - a3) .* ( (1-a3).*a3 );  
delta_2 = ( W2' * delta_3  + beta * repmat(sparsity_term, [1, m]) ).* ( (1-a2).*a2 );

W2grad = delta_3 * a2' ./ m + lambda * W2;
W1grad = delta_2 * a1' ./ m + lambda * W1;

b1grad = mean(delta_2,2);
b2grad = mean(delta_3,2);

sparse_cost = sparsityParam*log(sparsityParam./rho_hat) + (1-sparsityParam) * log((1-sparsityParam)./(1-rho_hat));
cost = sum(sum((a3-data).^2))/(2*m) + lambda/2 *(sum(sum(W1.^2))+sum(sum(W2.^2))) + beta * sum(sparse_cost);

%% -------------------------------------------------------------------

grad = [W1grad(:) ; W2grad(:) ; b1grad(:) ; b2grad(:)];

end

%-------------------------------------------------------------------

function sigm = sigmoid(x)
  
    sigm = 1 ./ (1 + exp(-x));
end

