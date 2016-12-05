load('generated_data/BRCA/Clinical.mat');
plotMSE('BRCA', 'Clinical', mse);
load('generated_data/BRCA/CNVGene.mat');
plotMSE('BRCA', 'CNVGene', mse);
load('generated_data/BRCA/Mutation.mat');
plotMSE('BRCA', 'Mutation', mse);
load('generated_data/BRCA/Protein.mat');
plotMSE('BRCA', 'Protein', mse);

load('generated_data/GBMLGG/Clinical.mat');
plotMSE('GBMLGG', 'Clinical', mse);
load('generated_data/GBMLGG/CNVGene.mat');
plotMSE('GBMLGG', 'CNVGene', mse);
load('generated_data/GBMLGG/Mutation.mat');
plotMSE('GBMLGG', 'Mutation', mse);
load('generated_data/GBMLGG/Protein.mat');
plotMSE('GBMLGG', 'Protein', mse);