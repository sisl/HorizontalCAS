# HorizontalCAS

This repository has three main components:
* Generate MDP table policy
* Compress MDP table policy by training a neural network
* Visualize MDP and neural network policies

## Generate MDP Policy
Required Julia Packages: Printf, POMDPs, POMDPModelTools, LocalFunctionApproximation, GridInterpolations, Distributed, SharedArrays, StaticArrays, HDF5

Tested with Julia v1.1

The policy is generated in parallel via Julia by running `julia -p NUM_PROCS SolveMatrix.jl` in the GenerateTable folder, where NUM_PROCS is the number of processors you want to use. The top of SolveMatrix.jl specifies where the table should be written to as an HDF5 file.

## Train Neural Networks
Required Python Packages: numpy, h5py, tensorflow 

Tested with Python v3.6

After generating the table, the table needs to be formatted into training data for the neural network. To do this, run `python genTrainingData.py` in the GenerateNetworks folder. The top of the file specifies the MDP table policy folder and the format for the resulting training data files. Note that the table is split into separate files, one for each previous advisory and tau combinations, which allows separate networks to be trained for each combination.

Next, run `python trainHCAS.py PREV_ADV TAU <gpu_ind>`, where PREV_ADV is the index of the previous advsiroy you want to train, TAU is the tau value, and gpu_ind is an optional input to specify which GPU to use. If you want to use a CPU instead, use -1 (the default if omitted). Options at the top of the file allow you to specify where the training data is stored and where the .nnet files should be written. There are additional options for the user to specify additional setting for network training.

## Visualize the Policies
Required Julia Packages: GridInterpolations, Interact, PGFPlots, Colors, ColorBrewer, HDF5

Tested with Julia v1.1

After generating MDP policies and training neural networks, the policies can be visualized. There is an example Jupyter notebook in the PolicyViz folder that shows how the policies can be interactively visualized.
