# Machine Learning Final Project

## How to setup environment:

### Clone this project to your local machine

```
$ git clone git@github.com:imgarylai/mlfp.git
```

### Setup data

1. Create data folder

```
cd mlfp
mkdir -p data/clean
```

2. Download and unzip [cleaned data](https://www.dropbox.com/sh/1njsvadxepv8dxa/AAB1FFMB_Fy9EOZFmTi5mSS-a?dl=1) to `mlfp/data/clean`

### Matlab setup

Open MATLAB and go to the navigation bar.

HOME -> ENVIRONMENT -> Set path ->Add with subdirectories

Add `mlfp` directory to MATLAB PATH

### File directory
```
mlfp/
|____ data/
|____ clean/
|     | cleanedData_BRCA_all.mat
|     | cleanedData_BRCA.mat
|____ example_code/
|     | cindex.m
|     | PerformanceExample.m
|____ extensions/
|     | ...
|____ GA/
|     | ...
|____ src/
|     |____ bin/
|           | ...
|     |____ feature_selections/
|           | ml_PCA.m
|           | test_feature_selection.m
|     |____ random_files/
|           | p-value test.m
|     |____ utils/
|           | cindex.m
|           | data_preprosssessing.m
|           | eva_cv.m
|           | eva.m
|     |____ demonstration.m

```
This setup guide helps you open our project in `mlfp` directory on MATLAB.

If you meet any problem, please let me know :) . (Gary)
