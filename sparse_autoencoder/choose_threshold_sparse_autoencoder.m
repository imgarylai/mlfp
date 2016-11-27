load data/clean/clean.mat;
mincost = 9999999;
for hiddenSize = 10:10:60
    for sparsityParam = 0:0.1:1
        for lambda = 0:0.1:1
            for beta = 1:1:10
                [newFeatures,cost] = ml_sparse_autoencoder(clinical_train.Features,hiddenSize,sparsityParam,lambda,beta);
                if cost < mincost
                    mincost = cost;
                    opt_hiddenSize = hiddenSize;
                    opt_beta = beta;
                    opt_lambda = lambda;
                    opt_sparsityParam = sparsityParam;
                end
            end
        end
    end
end
v.hiddenSize = opt_hiddenSize;
v.beta = opt_beta;
v.lambda = opt_lambda;
v.sparsityParam = opt_sparsityParam;
fname = sprintf('data/save/clinical_train');
save(fname, '-struct', 'v');
