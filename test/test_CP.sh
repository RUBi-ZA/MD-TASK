#!/bin/bash


mkdir out_CP
cd out_CP

BIN_DIR=../..

cp $BIN_DIR/example/* .

echo ""
echo "#### TEST COORDINATION PROPENSITY ####"
echo ""
python $BIN_DIR/cp.py --topology example_small.pdb --step 100 example_small.dcd
python $BIN_DIR/cp_analyse.py cp_cp_matrix.npy --avg --sliding_avg --window 4 --diff ${DIFF}
