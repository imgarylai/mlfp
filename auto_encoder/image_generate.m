load 'data/clean/clean.mat';
T = 5;
X = CNV_train.Features;
[featureSize, sampleSize] = size(X);
mse = zeros(T, featureSize);


for t = 1:T
    disp('-------------------------------------');
    fprintf('Interation: %d. \n', t);
    disp('-------------------------------------');
    for i=1:featureSize
        [newFeatures, autoenc, mseError]= ml_autoencoder(X, i);
        if(t == 1)
            [newFeatureSize, newSampleSize] = size(newFeatures);
            v.features = newFeatures;
            v.autoencoder = autoenc;
            fname = sprintf('auto_encoder/generated_data/autoencoder-%s-%d-%d.mat', char(CNV_train.SymbolTypes(1)) ,featureSize, newFeatureSize);
            save(fname, '-struct', 'v'); 
        end 
        mse(t, i) = mseError;
    end
end

meanMse = mean(mse);
stdMse = std(mse);
figure;
hold on;
title('MSE of the of the dimentional reduction by autoencoder');
xlabel('number of neuro in the hidden layer');
ylabel('Average MSE of 5 iterations');
H1 = plot(1:featureSize, meanMse, '--');
H2 = plot(1:featureSize, [meanMse - stdMse; meanMse + stdMse],':');
legend([H1,H2(1)],'mse','mse std', 'Location', 'Northeast');
hold off;