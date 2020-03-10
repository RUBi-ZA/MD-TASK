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


# Test output checksum
cd ..
check_sum=`find out_DCC -type f -exec md5sum {} \; | sort -k 2 | md5sum`

if [ "${check_sum}" != "e4e191a3989207cb5cd4c2847d069b1f  -" ]
then
    echo "Inconsistent output, this may be due to an error or modification to the algorithm or the inputs."
    echo "Please ensure that the output and checksum is correct."
    exit
fi
