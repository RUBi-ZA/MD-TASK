#!/bin/bash


mkdir out_PRS
cd out_PRS

BIN_DIR=../..

cp $BIN_DIR/example/* .

echo ""
echo "#### TEST PERTURBATION RESPONSE SCANNING ####"
echo ""
python $BIN_DIR/prs.py --final final.pdb --perturbations 20 --step 200 --prefix result --topology example_small.pdb example_small.xtc
