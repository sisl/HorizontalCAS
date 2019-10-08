## ReluVal Neural Network Over-approximation

This folder contains code for using ReluVal to over approximate the neural network policy.

To run the code (from this directory, on a Linux terminal with make installed), execute 
```bash
make 
mkdir Results
./runReluVal_HCAS.sh
```

The last command runs 10 instances of ReluVal in parallel at a time. 
Each ReluVal query is completely independent, so this process can be sped up if more than ten threads can be run at once.

The results are written to text files in the Results folder, which are read and used by the reachability analysis functions.
