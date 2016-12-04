clear
addpath data/
addpath CNV/
addpath clinical/
addpath mutation/
addpath protein/
addpath mRNA/

warning('off','all')
original_BRCA = load ('GBMLGG.Data.mat');
prepro_BRCA=prepro(original_BRCA);

Clinical= getAvailableClinical(prepro_BRCA);
Clinical= rmirrelevant(Clinical);

CNV= getAvailableCNV(prepro_BRCA);
CNV= rmirrelevant(CNV);

Mutation= getAvailableMutation(prepro_BRCA);
Mutation= rmirrelevant(Mutation);

Protein= getAvailableProtein(prepro_BRCA);
Protein= rmirrelevant(Protein);


