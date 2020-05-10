# HorizontalCAS

This repository describes how to generate HorizontalCAS score tables and train a neural network representation. HorizontalCAS is a simple, notional collision avoidance system (CAS) that gives horizontal turning advisories to an aircraft to avoid an intruder. This system is inspired by early prototypes of ACAS Xu and neural networks trained to represent the score table, but HorizontalCAS is not related in any way to ACAS. This simple system is inteded to facilitate research towards safety-critical neural networks and their verification in the hopes that lessons learned can be applied to real systems.

This repository supports a paper presented at the Digital Avionics Systems Conference (DASC) in 2019, which can be found [here](https://arxiv.org/pdf/1912.07084.pdf).

If HorizontalCAS is useful for your research, please cite
```
@inproceedings{julian2019guaranteeing,
  title={Guaranteeing safety for neural network-based aircraft collision avoidance systems},
  author={Julian, Kyle D and Kochenderfer, Mykel J},
  booktitle={Digital Avionics Systems Conference (DASC)},
  year={2019}
}
```

## Saved training data and neural networks
Although this repository contains the source code used to generate the advisory score table and train a neural network representation, saved copies of training data and trained neural networks can be found here: 
* [**Training Data**](https://drive.google.com/drive/folders/14kcGM_G5sq72BpCfD4dimp27S7ael3by?usp=sharing)
* [**Trained Networks**](https://drive.google.com/drive/folders/1Sj2noNh65xbG6H1fO3DkS1GnevSYTa5b?usp=sharing)

The remainder of this README describes how the score table is generated in Julia, how neural networks are trained in Python using tensorflow, and how the neural network policies can be visualized using Julia kernel for a Jupyter notebook.

## Generate MDP Policy
Required Julia Packages: Printf, POMDPs, POMDPModelTools, LocalFunctionApproximation, GridInterpolations, Distributed, SharedArrays, StaticArrays, HDF5

Tested with Julia v1.1

The policy is generated in parallel via Julia by running `julia -p NUM_PROCS SolveMDP.jl` in the GenerateTable folder, where NUM_PROCS is the number of processors you want to use. The top of SolveMatrix.jl specifies where the table should be written to as an HDF5 file.

## Train Neural Networks
Required Python Packages: numpy, h5py, tensorflow 

Tested with Python v3.6

After generating the table, the table needs to be formatted into training data for the neural network. To do this, run `python genTrainingData.py` in the GenerateNetworks folder. The top of the file specifies the MDP table policy folder and the format for the resulting training data files. Note that the table is split into separate files, one for each previous advisory and tau combinations, which allows separate networks to be trained for each combination.

Next, run `python trainHCAS.py PREV_ADV TAU <gpu_ind>`, where PREV_ADV is the index of the previous advsiroy you want to train, TAU is the tau value, and gpu_ind is an optional input to specify which GPU to use. If you want to use a CPU instead, use -1 (the default if omitted). Options at the top of the file allow you to specify where the training data is stored and where the .nnet files should be written. There are additional options for the user to specify additional setting for network training.

## Visualize the Policies
Required Julia Packages: GridInterpolations, Interact, PGFPlots, Colors, ColorBrewer, HDF5

Tested with Julia v1.1

After generating MDP policies and training neural networks, the policies can be visualized. There is an example Jupyter notebook in the PolicyViz folder that shows how the policies can be interactively visualized.
