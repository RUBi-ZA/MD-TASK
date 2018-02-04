#!/bin/bash


mkdir out_L
cd out_L

BIN_DIR=../..

cp $BIN_DIR/example/* .

echo ""
echo "#### AVERAGE SHORTEST PATH - WT ####"
echo ""

PREFIX=wt

python $BIN_DIR/calc_network.py --topology $PREFIX.pdb --threshold 7.0 --step 100 --generate-plots --calc-L --discard-graphs $PREFIX.dcd

mv ${PREFIX}_0_avg_L.dat ref_${PREFIX}_avg_L.dat

python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat --normalize
python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat --normalize --normalization-mode nonzero
python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat --normalize --normalization-mode plusone
python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat --normalize --normalization-mode standard
python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat
python $BIN_DIR/avg_network.py --data ${PREFIX}_*_avg_L_delta_L.dat --data-type delta-L --prefix ${PREFIX} --generate-plots --x-label "Residues" --y-label "Avg delta L" --title "Wild Type"


echo ""
echo ""
echo "#### AVERAGE SHORTEST PATH - MUTANT ####"
echo ""

PREFIX=mutant

python $BIN_DIR/calc_network.py --topology $PREFIX.pdb --threshold 7.0 --step 100 --generate-plots --calc-L --discard-graphs --lazy-load $PREFIX.dcd

mv ${PREFIX}_0_avg_L.dat ref_${PREFIX}_avg_L.dat

python $BIN_DIR/calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat --normalize
python $BIN_DIR/avg_network.py --data ${PREFIX}_*_avg_L_delta_L.dat --data-type delta-L --prefix ${PREFIX} --generate-plots --x-label "Residues" --y-label "Avg delta L" --title "Mutant"


echo ""
echo ""
echo "#### AVERAGE SHORTEST PATH - WT vs MUTANT ####"
echo ""

python $BIN_DIR/compare_networks.py --prefix "wt_mutant_L_avg" --reference-label Wild-type --alternative-label Mutant --y-label "Delta L" --reference wt_delta_L_avg.dat --alternative mutant_delta_L_avg.dat
python $BIN_DIR/compare_networks.py --prefix "wt_mutant_L_std_dev" --reference-label Wild-type --alternative-label Mutant --y-label "Delta L" --reference wt_delta_L_std_dev.dat --alternative mutant_delta_L_std_dev.dat
python $BIN_DIR/delta_networks.py --reference wt_delta_L_avg.dat --reference-std wt_delta_L_std_dev.dat --alternatives *delta_L_avg.dat --alternatives-std *delta_L_std_dev.dat --absolute --prefix wt_mutant_cmap --title "My Protein" --x-label "Residues" --y-label "Proteins"
