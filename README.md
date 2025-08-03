This repository contains source code for fitting a Gaussian DPMM using a combination of algorithms described in Neal (2000) and Jain & Neal (2004, 2007).

Currently, five different model structures are supported:
- Conjugate DEE, DEV, and UVV
- Non-conjugate DEV, UVV

Please see the example files at the top level of this repo.

The models require the following inputs:
- Data y: this is a length n list of p-dimensional vectors for each observation in the data set on which the user wishes to perform clustering
- Prior hyperparameters as needed for the specified type of model
- An initial number of groups, k_init
