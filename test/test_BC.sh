#!/bin/bash


mkdir out_BC
cd out_BC

BIN_DIR=../..

cp $BIN_DIR/example/* .

echo ""
echo "#### BETWEENNESS CENTRALITY - WT ####"
echo ""

PREFIX=wt
NORM=plusone

python $BIN_DIR/calc_network.py --topology $PREFIX.pdb --threshold 7.0 --step 100 --generate-plots --calc-BC --save-graphs --lazy-load $PREFIX.dcd --calc-BC

mv ${PREFIX}_0_BC.dat ref_${PREFIX}_BC.dat

python $BIN_DIR/calc_delta_BC.py --generate-plots --normalize --normalization-mode ${NORM} --reference ref_${PREFIX}_bc.dat --alternatives ${PREFIX}_*_bc.dat
python $BIN_DIR/avg_network.py --data ${PREFIX}_*_bc_${NORM}_delta_BC.dat --data-type delta-BC --prefix ${PREFIX} --generate-plots --x-label "Residues" --y-label "Avg delta BC" --title "Wild Type"


echo ""
echo ""
echo "#### BETWEENNESS CENTRALITY - MUTANT ####"
echo ""

PREFIX=mutant

python $BIN_DIR/calc_network.py --topology $PREFIX.pdb --threshold 7.0 --step 100 --generate-plots --calc-BC --save-graphs $PREFIX.dcd --calc-BC

mv ${PREFIX}_0_BC.dat ref_${PREFIX}_BC.dat

python $BIN_DIR/calc_delta_BC.py --generate-plots --normalize --normalization-mode ${NORM} --reference ref_${PREFIX}_bc.dat --alternatives ${PREFIX}_*_bc.dat
python $BIN_DIR/avg_network.py --data ${PREFIX}_*_bc_${NORM}_delta_BC.dat --data-type delta-BC --prefix ${PREFIX} --generate-plots --x-label "Residues" --y-label "Avg delta BC" --title "Mutant"


echo ""
echo ""
echo "#### BETWEENNESS CENTRALITY - WT vs MUTANT ####"
echo ""

python $BIN_DIR/compare_networks.py --prefix "wt_mutant_BC_avg" --reference-label Wild-type --alternative-label Mutant --y-label "Delta BC" --reference wt_delta_BC_avg.dat --alternative mutant_delta_BC_avg.dat
python $BIN_DIR/compare_networks.py --prefix "wt_mutant_std_dev" --reference-label Wild-type --alternative-label Mutant --y-label "Delta BC" --reference wt_delta_BC_std_dev.dat --alternative mutant_delta_BC_std_dev.dat
python $BIN_DIR/delta_networks.py --reference wt_delta_BC_avg.dat --reference-std wt_delta_BC_std_dev.dat --alternatives *delta_BC_avg.dat --alternatives-std *delta_BC_std_dev.dat --absolute --prefix wt_mutant_cmap --title "My Protein" --x-label "Residues" --y-label "Proteins"


