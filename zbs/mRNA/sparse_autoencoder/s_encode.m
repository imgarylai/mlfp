function answer = s_encode (weight,b ,input)

W1 = weight;
b1 = b;
z2 = W1 * input + repmat(b1,[1, size(input,2)]);
answer = sigmoid(z2);