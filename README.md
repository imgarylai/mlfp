# Machine Learning Final Project

## How to setup environment:

### Clone this project to your local machine

```
$ git clone git@github.com:imgarylai/mlfp.git
```

### Setup data

Download and unzip [cleaned data](https://github.com/imgarylai/mlfp/wiki/Clean-data-%5BFinal%5D) to `mlfp/data/`

### Matlab setup

Open MATLAB and go to the navigation bar.

HOME -> ENVIRONMENT -> Set path ->Add with subdirectories

Add `mlfp` directory to MATLAB PATH

### How to run

Run main.m, then you will be prompted to enter. 

### File directory
```
mlfp/
|____ main.m
|
|____ cIndex.m
|
|____ data/
|     | BRCA.mat (download needed)
|     | GBMLGG.mat (download needed)
|
|____ autoencoder/
|     |autoencoder.m
|     
|____ genetic_algorithm/
|     | c_index_fitness.m
|     | PopFunction.m
|     
|____ PCA/
|     | newPCA.m
|     |____ extensions/
|           | ...
|
|____ preprocessing/
|     | ...
|
|____ sparse_autoencoder/
|     | newPCA.m
|     | ...
|     |____ minFunc/
|           | ...
|
```
This setup guide helps you open our project in `mlfp` directory on MATLAB.

If you meet any problem, please let me know :) . (Gary)
