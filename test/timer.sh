#!/bin/bash

PREFIX=example_small

mkdir timer
cd timer

cp ../../example/* .


echo ""
echo ""
echo "#### TIMING: calc_network.py (--calc-BC)"
echo ""

for i in `seq 1 10`;
do
    time python ../../calc_network.py --topology ${PREFIX}.pdb --threshold 7.0 --step 100 --generate-plots --calc-BC --discard-graphs ${PREFIX}.dcd > /dev/null
done  



echo ""
echo ""
echo "#### TIMING: calc_network.py (--calc-L)"
echo ""

for i in `seq 1 10`;
do
    time python ../../calc_network.py --topology example_small.pdb --threshold 7.0 --step 100 --generate-plots --calc-L --discard-graphs example_small.dcd > /dev/null
done  



echo ""
echo ""
echo "#### TIMING: calc_delta_BC.py"
echo ""

mv ${PREFIX}_0_bc.dat ref_${PREFIX}_bc.dat 

for i in `seq 1 10`;
do
    time python ../../calc_delta_BC.py --generate-plots --reference ref_${PREFIX}_bc.dat --alternatives ${PREFIX}_*_bc.dat > /dev/null
done   



echo ""
echo ""
echo "#### TIMING: calc_delta_L.py"
echo ""

mv ${PREFIX}_0_avg_L.dat ref_${PREFIX}_avg_L.dat 

for i in `seq 1 10`;
do
    time python ../../calc_delta_L.py --generate-plots --reference ref_${PREFIX}_avg_L.dat --alternatives ${PREFIX}_*_avg_L.dat > /dev/null
done 



echo ""
echo ""
echo "#### TIMING: avg_network.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../avg_network.py --data ${PREFIX}_*_bc_delta_BC.dat --data-type delta-BC --prefix ${PREFIX} --generate-plots --x-label "Residues" --y-label "Avg delta BC" --title "Wild Type" > /dev/null
done 



echo ""
echo ""
echo "#### TIMING: compare_networks.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../compare_networks.py --prefix "${PREFIX}_BC_avg" --reference-label Wild-type --alternative-label Wild-type-2 --y-label "Delta BC" --reference ${PREFIX}_delta_BC_avg.dat --alternative ${PREFIX}_delta_BC_avg.dat > /dev/null
done 



echo ""
echo ""
echo "#### TIMING: delta_networks.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../delta_networks.py --reference ${PREFIX}_delta_BC_avg.dat --reference-std ${PREFIX}_delta_BC_std_dev.dat --alternatives ${PREFIX}_delta_BC_avg.dat ${PREFIX}_delta_BC_avg.dat --alternatives-std ${PREFIX}_delta_BC_std_dev.dat ${PREFIX}_delta_BC_std_dev.dat --absolute --prefix ${PREFIX}_cmap --title "My Protein" --x-label "Residues" --y-label "Proteins" > /dev/null
done 



echo ""
echo ""
echo "#### TIMING: contact_map.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../contact_map.py --residue LEU9 --prefix ${PREFIX} --topology ${PREFIX}.pdb ${PREFIX}.dcd > /dev/null
done  



echo ""
echo ""
echo "#### TIMING: calc_correlation.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../calc_correlation.py --step 100 --prefix example_corr --trajectory ${PREFIX}.dcd --topology ${PREFIX}.pdb --lazy-load > /dev/null
done   



echo ""
echo ""
echo "#### TIMING: prs.py"
echo ""

for i in `seq 1 10`;
do
    time python ../../prs.py --initial initial.xyz --final final.xyz --perturbations 50 --step 100 --prefix result --topology ${PREFIX}.pdb ${PREFIX}.dcd > /dev/null
done  






