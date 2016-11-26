opttheta = [];
optbeta = 0;
mincost = 9999999;
for hiddenSize = 10:10:60
    for sparsityParam = 0:0.1:1
        for lambda = 0:0.1:1
            for beta = 1:1:10
                [newFeatures,cost] = ml_sparse_autoencoder(clinical_train.Features,hiddenSize,sparsityParam,lambda,beta);
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
end
    

