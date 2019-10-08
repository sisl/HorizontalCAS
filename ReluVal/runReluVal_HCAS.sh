#!/bin/bash

# There are 5 previous RA (pra) and 8 tau combinations for 40 networks total
# This script runs ReluVal to compute an over-approximation of batches of ten networks at a time
# The results are saved in a text file, which can be read by functions in the Reachability directory

for pra in `seq 0 4`;
do
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau00_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau00.txt &
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau05_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau05.txt &
done
wait

for pra in `seq 0 4`;
do
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau10_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau10.txt &
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau15_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau15.txt &
done
wait

for pra in `seq 0 4`;
do
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau20_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau20.txt &
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau30_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau30.txt &
done
wait

for pra in `seq 0 4`;
do
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau40_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau40.txt &
    ./network_test_HCAS ../networks/HCAS_rect_v6_pra${pra}_tau60_25HU_3000.nnet > ./Results/HCAS_v6_pra${pra}_tau60.txt &
done
wait

