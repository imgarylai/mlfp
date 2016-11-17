tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% load data
clc; clear;
load '../data/BRCA.Data.mat';

maxNumCompThreads(8);
% warning('off','all')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% transform data
Symbols         =   cellstr(Symbols)';
SymbolTypes     =   cellstr(SymbolTypes)';
Survival        =   Survival';
Censored        =   Censored';
Features        =   Features';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% separate features according to their types, part 1
clinical_indices    = strcmp(SymbolTypes, 'Clinical');
CNV_indices         = strcmp(SymbolTypes, 'CNVGene') | strcmp(SymbolTypes, 'CNVArm');
mutation_indices    = strcmp(SymbolTypes, 'Mutation');
mRNA_indices        = strcmp(SymbolTypes, 'mRNA');
protein_indices     = strcmp(SymbolTypes, 'Protein');
assert(sum(clinical_indices) + sum(CNV_indices) + sum(mutation_indices) ...
       + sum(mRNA_indices) + sum(protein_indices) == size(Features, 2), ...
       'after split features, sum not equal');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drop samples with Nan in Censored or Survival
fprintf('number of Nan in Censored %d\n', sum(isnan(Censored)));
fprintf('number of Nan in Survival %d\n', sum(isnan(Survival)));
sample_keep_indices = ~(isnan(Censored) | isnan(Survival));
Features    =   Features(sample_keep_indices, :);
Survival    =   Survival(sample_keep_indices);
Censored    =   Censored(sample_keep_indices);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drop samples with only Nan in mRNA and protein
features_mRNA       = Features(:, mRNA_indices);
features_protein    = Features(:, protein_indices);
features_mRNA_protein = [features_mRNA, features_protein];

fprintf('number of samples with only Nan in mRNA and protein %d\n', ...
        sum(all(isnan(features_mRNA_protein), 2)));
sample_keep_indices = ~all(isnan(features_mRNA_protein), 2);

Features    =   Features(sample_keep_indices, :);
Survival    =   Survival(sample_keep_indices);
Censored    =   Censored(sample_keep_indices);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% separate features according to their types, part 2, because some
%% samples are dropped, so this needs to be done after dropping samples
features_clinical   = Features(:, clinical_indices);
features_CNV        = Features(:, CNV_indices);
features_mutation   = Features(:, mutation_indices);
features_mRNA       = Features(:, mRNA_indices);
features_protein    = Features(:, protein_indices);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% find samples with Nan in Clinical or CNV or Mutation

% fprintf('number of samples missing Clinical %d\n', sum(strcmp(AvailableClinical, 'No')))
% fprintf('number of samples missing CNV %d\n', sum(strcmp(AvailableCNV, 'No')))
% fprintf('number of samples missing Mutation %d\n', sum(strcmp(AvailableMutation, 'No')))

fprintf('number of samples missing Clinical %d\n', ...
        sum( (sum(isnan(Features(:, clinical_indices)), 2)==0)==0))
fprintf('number of samples missing CNV %d\n', ...
        sum( (sum(isnan(Features(:, CNV_indices)), 2)==0)==0))
fprintf('number of samples missing Mutation %d\n', ...
        sum( (sum(isnan(Features(:, mutation_indices)), 2)==0)==0))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% make new features indicating whether given feature is Nan or not

% generate new features/indicators
indicator_clinical  = isnan(features_clinical);
indicator_CNV       = isnan(features_CNV);
indicator_mutation  = isnan(features_mutation);
% clear Nan in original features to 0
features_clinical(indicator_clinical)   = 0;
features_CNV(indicator_CNV)             = 0;
features_mutation(indicator_mutation)   = 0;
% combine cleared features and indicator features
features_clinical   =   [features_clinical, ~indicator_clinical];
features_CNV        =   [features_CNV, ~indicator_CNV];
features_mutation   =   [features_mutation, ~indicator_mutation];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use knnimputation for Nan in mRNA and protein

DROP_FEATURE_THRESHOLD  = 1/3;
IMPUTATION_K = 3;
features_mRNA_protein = [features_mRNA, features_protein];

% step 1. within a feature, if number of missing(Nan) larger than
% total records * DROP_FEATURE_THRESHOLD, then drop current feature
fprintf('before dropping features with too many Nan, number of features %d\n', size(features_mRNA_protein, 2))
feature_keep_indices    = (sum(isnan(features_mRNA_protein), 1)) < ...
length(Censored) * DROP_FEATURE_THRESHOLD;
features_mRNA_protein   = features_mRNA_protein(:, feature_keep_indices);
% Symbols = Symbols(feature_keep_indices);
% SymbolTypes = SymbolTypes(feature_keep_indices);
fprintf('after dropping features with too many Nan, number of features %d\n\n', size(features_mRNA_protein, 2))

% step 2. now drop all features with range == 0
fprintf('before dropping features with range 0, number of features %d\n', size(features_mRNA_protein, 2))
feature_keep_indices    = (range(features_mRNA_protein, 1) ~= 0);
features_mRNA_protein   = features_mRNA_protein(:, feature_keep_indices);
% Symbols = Symbols(feature_keep_indices);
% SymbolTypes = SymbolTypes(feature_keep_indices);
fprintf('after dropping features with range 0, number of features %d\n\n', size(features_mRNA_protein, 2))

% step 3. imputation using euclidean distance with k neighbours
fprintf('before imputation, number of Nan %d\n', sum(sum(isnan(features_mRNA_protein))));
% features_mRNA_protein = knnimpute(features_mRNA_protein, IMPUTATION_K, 'Distance', 'seuclidean');
features_mRNA_protein(isnan(features_mRNA_protein)) = 0;
fprintf('after imputation, number of Nan %d\n\n', sum(sum(isnan(features_mRNA_protein))));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% normalize mRNA and protein features(zscore) and age in clinical
% features_mRNA_protein_cov   = repmat(diag(cov(features_mRNA_protein))', ...
%                                  size(features_mRNA_protein, 1), 1);
% features_mRNA_protein_cov_mean = repmat(mean(features_mRNA_protein, 1), ...
%                                     size(features_mRNA_protein, 1), 1);
% features_mRNA_protein = (features_mRNA_protein - features_mRNA_protein_cov_mean);
% features_mRNA_protein = features_mRNA_protein./features_mRNA_protein_cov;
features_mRNA_protein = zscore(features_mRNA_protein);
features_clinical(:, 1) = zscore(features_clinical(:, 1));

toc;
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% save results
save('../data/cleanedData_BRCA.mat', 'features_clinical', 'features_CNV', 'features_mutation', ...
     'features_mRNA_protein', 'Survival', 'Censored')
save('../data/cleanedData_BRCA_all.mat')

toc;
