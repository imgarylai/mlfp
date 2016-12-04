function image_generator(prefix ,type)
    T = 5;
    X = type.Features;
    [featureSize, ~] = size(X);
    mse = zeros(T, featureSize);
    save_path = sprintf('generated_data/%s/%s.mat', char(prefix), char(type.SymbolTypes(1)));

    for t = 1:T
        disp('-------------------------------------');
        fprintf('Interation: %d. \n', t);
        disp('-------------------------------------');
        for i=1:featureSize
            [newFeatures, autoenc, mseError]= ml_autoencoder(X, i);
    %         if(t == 1)
    %             [newFeatureSize, newSampleSize] = size(newFeatures);
    %             v.features = newFeatures;
    %             v.autoencoder = autoenc;
    %             fname = sprintf('generated_data/autoencoder-%s-%d-%d.mat', char(type.SymbolTypes(1)) ,featureSize, newFeatureSize);
    %             save(fname, '-struct', 'v'); 
    %         end 
            mse(t, i) = mseError;
            save(save_path , 'mse');
        end
    end

end

