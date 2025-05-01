# Tensor Decomposition Projects

Data projects based on the book *Tensor Decompositions for Data Science* by Ballard and Kolda (March 2025 edition).

This repository contains implementations and comparisons of various tensor decomposition methods.

## Methods Implemented

- **Tucker Decomposition**
  - HOSVD
  - ST-HOSVD
  - HOOI
- **CP Decomposition**
  - CP-ALS
  - CP-OPT
  - CP-NLS

## Folder Structure

- `code/`: MATLAB implementations of the algorithms
- `figures/`: Plots and visualizations of decomposition results
- `data/`: Sample tensors (synthetic and real-world)
- `README.md`: Project overview

## Datasets

The `data/` folder contains sample tensors used to test and compare decomposition algorithms. These include:

- **Synthetic tensors** generated from known factor matrices, with and without added noise.
- **Real-world example**: A tensor representing Chicago crime data  
  *(dimensions: days × hours × community areas × crime types)*  
  used to test compression and anomaly detection performance.

If you use your own datasets, please make sure they follow the same format (`.mat` or `.csv`) and adjust the scripts accordingly.

## Getting Started

To run the code:
1. Open the `.m` files in MATLAB
2. Load data from the `data/` folder as instructed in each script
3. View results and plots in the `figures/` folder

## Requirements

All code in this repository requires the **Tensor Toolbox for MATLAB** by Bader and Kolda.  
You can download it from the official website:

🔗 [https://www.tensortoolbox.org](https://www.tensortoolbox.org)

After downloading:
1. Unzip the folder
2. Add it to your MATLAB path