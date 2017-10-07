#!/bin/bash

PREFIX=example_small

mkdir -p timer/prs
cd timer/prs

cp ../../../example/* .

echo ""
echo ""
echo "#### TIMING: prs.py"
echo ""

for i in $(seq 1 10);
do
    time python ../../../prs.py --initial initial.xyz --final final.xyz --perturbations 50 --step 100 --prefix result --topology ${PREFIX}.pdb ${PREFIX}.dcd > /dev/null
done
