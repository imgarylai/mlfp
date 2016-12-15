clear
warning('off','all')
addpath ..
BRCA_Original = load ('../data/BRCA.Data.mat');
GBMLGG_Original = load ('../data/GBMLGG.Data.mat');

BRCA_prepro = prepro(BRCA_Original);
GBMLGG_prepro = prepro(GBMLGG_Original);

% BRCA.Clinical=getAvailableClinical(BRCA_prepro);
% BRCA.Mutation=getAvailableMutation(BRCA_prepro);
% BRCA.CNV=getAvailableCNV(BRCA_prepro);
% BRCA.Protein=getAvailableProtein(BRCA_prepro);
BRCA.mRNA=getAvailablemRNA(BRCA_prepro);

% GBMLGG.Clinical=getAvailableClinical(GBMLGG_prepro);
% GBMLGG.Mutation=getAvailableMutation(GBMLGG_prepro);
% GBMLGG.CNV=getAvailableCNV(GBMLGG_prepro);
% GBMLGG.Protein=getAvailableProtein(GBMLGG_prepro);
GBMLGG.mRNA=getAvailablemRNA(GBMLGG_prepro);


% BRCA.Clinical=rmirrelevant(BRCA.Clinical);
% BRCA.Mutation=rmirrelevant(BRCA.Mutation);
% BRCA.CNV=rmirrelevant(BRCA.CNV);
% BRCA.Protein=rmirrelevant(BRCA.Protein);
BRCA.mRNA=rmirrelevant(BRCA.mRNA);

% GBMLGG.Clinical=rmirrelevant(GBMLGG.Clinical);
% GBMLGG.Mutation=rmirrelevant(GBMLGG.Mutation);
% GBMLGG.CNV=rmirrelevant(GBMLGG.CNV);
% GBMLGG.Protein=rmirrelevant(GBMLGG.Protein);
GBMLGG.mRNA=rmirrelevant(GBMLGG.mRNA);