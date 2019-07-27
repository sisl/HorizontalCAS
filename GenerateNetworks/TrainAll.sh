#!/bin/bash
for pra in 0 1 2 3 4
do
    python trainHCAS.py $pra 0 0 &
    python trainHCAS.py $pra 5 0 &
    python trainHCAS.py $pra 10 1 &
    python trainHCAS.py $pra 15 1 &
    python trainHCAS.py $pra 20 2 &
    python trainHCAS.py $pra 30 2 &
    python trainHCAS.py $pra 40 3 &
    python trainHCAS.py $pra 60 3 &
    wait
done
