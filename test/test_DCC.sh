#!/bin/bash


mkdir out_DCC
cd out_DCC

BIN_DIR=../..

cp $BIN_DIR/example/* .

echo ""
echo "#### TEST DYNAMIC CROSS CORRELATION ####"
echo ""

python $BIN_DIR/calc_correlation.py --step 100 --prefix example_corr --trajectory wt.dcd --topology wt.pdb --lazy-load
python $BIN_DIR/calc_correlation.py --step 100 --prefix example_corr --trajectory mutant.dcd --topology mutant.pdb --lazy-load
