# Data exploration in the sRNA datasets

## Environment setup

```bash
cd THIS_PATH_HERE

# setup
conda create python=3.10 --prefix ./.condaenv -y

# start the environment
conda activate ./.condaenv

# dependencies
conda install -y matplotlib pandas seaborn numpy biopython
conda install -y jupyter
```

## Execution

To run this project, follow these steps:

1. You need to obtain the datasets by running the `download_*.ipynb` files;
2. Feature extraction (`feature_extraction.ipynb`);
3. Data exploration (`data_exploration.ipynb`).
