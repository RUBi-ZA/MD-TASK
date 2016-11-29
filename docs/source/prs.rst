Pertubation Response Scanning
===============================

PRS is a computational technique that is useful for determining single residues that play an active role in the manipulation of protein conformational changes. As input it requires two distinct atomic conformations for a protein of interest; initial and target structures respectively. The technique then performs a residue-by-residue scanning of the initial conformation, by exerting multiple factious external forces of both random direction and magnitude on each residue in the protein structure. After external force perturbation, the subset of residues/forces that invoke a conformational change closest to the target structure are recorded. To calculate the predicted displacement of all residues in relation to a perturbation at a single residue, PRS requires the construction of a variance-covariance matrix, which can be obtained from suitable length MD simulation trajectories of the initial protein structure. The quality of the predicted displacements is then assessed by correlating the predicted and experimental displacements, averaged over all affected residues. This results in a correlation coefficient for each residue in the protein, where a value close to 1 implies good agreement with the experimental change. PRS can thus be used to map regions on a protein whose perturbation leads to a conformational change that resembles the expected target structure. These regions are often active site residues on the protein, but also potentially point to locations involved in allostery and allosteric control. PRS has also been used in conjunction with molecular docking to calculate ligand bound conformations from an unbound structure, in a scheme for exploring protein-ligand interaction.

Performing PRS
---------------

**Command:** :: 
	
	python prs.py <options> --final <final.xyz> --trajectory <trajectory> --topology <pdb>

**Inputs:**

===========================  ===========  ====================  ===========================================================================================================================================================================
 Input (*\*required*)        Input type   Flag                  Description                  
===========================  ===========  ====================  ===========================================================================================================================================================================
Trajectory *                 File                               A trajectory from a molecular dynamics simulation. Can be in DCD or XTC format.
Topology *                   File         ``--topology``        A PDB reference file for the trajectory.
Initial                      File         ``--initial``         Co-ordinate file (.xyz) depicting the initial conformation (default: co-ordinate file is generated from the first frame of the trajectory)
Final *                      File         ``--final``           Co-ordinate file (.xyz) depicting the target conformation
Perturbations                Integer      ``--perturbations``   Number of perturbations to apply
No. of frames in trajectory  Integer      ``--num-frames``      Optionally specify the number of frames in the trajectory. This will run the script in a memory efficient mode. Usefult for large trajectories that don't fit into memory.
Step                         Integer      ``--step``            Step to use when iterating through trajectory frames i.e. how many frames will be skipped.
Prefix                       Text         ``--prefix``          Prefix used to name outputs 
===========================  ===========  ====================  ===========================================================================================================================================================================

Given a trajectory, ``example_small.dcd``, with initial and target co-odinate files, ``initial.xyz`` and ``final.xyz``, respectively, and topology file, ``example_small.pdb``, the following command could be used: ::

	python prs.py --initial initial.xyz --final final.xzy --perturbations 250 --step 100 --prefix result --trajectory example_small.dcd --topology example_small.pdb


**Outputs:**

=====================  ===================================================================================================================================================================
Output                 Description
=====================  ===================================================================================================================================================================
Correlation CSV file   Correlation coefficient for each residue in the protein, where a value close to 1 implies good agreement with the experimental change
=====================  ===================================================================================================================================================================
