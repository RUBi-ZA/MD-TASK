Dynamic Cross-Correlation
=========================

Molecular Dynamics (MD) is a computational method that analyses the physical motions of atoms within a protein or protein complex. In a given system, the interactions between the atoms can be simulated in the presence of a force field and, following the application of Newtons' equations of motion, trajectories corresponding to the dynamical motions of the atoms are obtained. The trajectories represent sequential snapshots of the system, by presenting the atomic coordinates at specific time intervals throughout the simulation. This allows for the investigation into the dynamical changes of the system over time. The applications of MD simulations are vast. By analysing the trajectories of the system, it is possible to calculate the dynamic correlation between all atoms within the molecule i.e. the degree to which they move together. This dynamic cross-correlation tool produces an NxN heatmap, where N = the number of (alpha carbon) atoms in the system, and each element corresponds to the dynamic cross-correlation between each i,j atom. The correlation values are calculated between -1 and 1, where 1=complete correlation; -1=complete anti-correlation; 0= no correlation. 

Calculating dynamic cross-correlation
---------------------------------------

**Command:** :: 
	
	calc_correlation.py <options> --trajectory <trajectory> --topology <pdb>

**Inputs:**

=========================  ===========  ====================  ========================================================================================================================================================
 Input (*\*required*)      Input type   Flag                  Description                  
=========================  ===========  ====================  ========================================================================================================================================================
Trajectory *               File                               A trajectory from a molecular dynamics simulation. Can be in DCD or XTC format.
Topology *                 File         ``--topology``        A PDB reference file for the trajectory.
Step                       Integer      ``--step``            Step to use when iterating through trajectory frames i.e. how many frames will be skipped.
Prefix                     Text         ``--prefix``          Prefix used to name outputs.
Lazy load                  Boolean      ``--lazy-load``       Load trajectory frames in a memory efficient manner - use for large trajectories.
=========================  ===========  ====================  ========================================================================================================================================================

Given a trajectory, ``example_small.dcd``, and topology file, ``example_small.pdb``, the following command could be used: ::

	calc_correlation.py --step 100 --prefix example_corr --trajectory example_small.dcd --topology example_small.pdb --lazy-load



**Outputs:**

=====================  ===================================================================================================================================================================
Output                 Description
=====================  ===================================================================================================================================================================
Correlation heatmap    PNG heatmap depicting the dynamic correlation between atoms in the trajectory 
Correlation text file  Correlation data in text format
=====================  ===================================================================================================================================================================
