Network Analysis
==================

Residue Interaction Networks (RIN) are analyzed using a branch of Mathematics known as graph theory. In a RIN, each residue in the protein is a node in the network. An edge (or connection) between two nodes exists if there is an interaction between the two residues those nodes represent. MD-RIN considers an interaction between two residues to exist if the beta carbon atoms of the residues are within a user-defined cut-off distance (usually around 6.5 – 7.5 Å) of each other. Once the network has been constructed, there are various network measures that can be used to analyze it. Currently, MD-RIN can be used to analyze the change in betweenness centrality (BC) and average shortest path (L) of residues in a protein over a molecular dynamics simulation. This can be used to determine which residues are important for intra-protein communication and conformational changes. RINs can also be useful in the analysis of SNPs. Comparing changes in BC and L between the simulation of a wild-type and mutant protein can provide interesting insights into differences in intra-protein communication, which can affect the function of the protein.

Measurements
-----------------

**1. Betweenness Centrality (BC)**

Betweenness centrality (BC) is a measure of how important a residue is for communication within a protein. It is equal to the number of shortest paths from all vertices to all others that pass through that node. Residues in a protein that have a high BC reveal locations that tend to be important for controlling inter-domain commuication in a protein.

**2. Average Shortest Path(L)**

The average shortest path (L) to a given residue is calculated by working out all the shortest paths to the given node and dividing by the number of paths. The average shortest path to a residue gives an idea of how accessible the residue is within the protein. This can be used to, for example, analyze SNPs. A mutation may result in a change in L of a number of residues in the protein. This may indicate that the mutation has an important effect on protein function e.g. previous studies have suggested that positions that result in high delta L values may steer conformational changes.

Calculating BC and L
------------------------

**Command:** ::
	
	python calc_network.py <options> --topology <pdb file> <trajectory>

**Inputs:**

+------------------------+------------+--------------------+-----------------------------+
| Input (*\*required*)   | Input type | Flag               | Description                 |
+========================+============+====================+=============================+
|Trajectory *            | File       |                    |A trajectory from a molecular| 
|                        |            |                    |dynamics simulation. Can be  |
|                        |            |                    |in DCD or XTC format.        |
+------------------------+------------+--------------------+-----------------------------+
|Topology *              | File       |``--topology``      |A PDB reference file for the |
|                        |            |                    |trajectory.                  |
+------------------------+------------+--------------------+-----------------------------+ 
|Ligands                 | CSV ligand |``--ligands``       |Ligands that should be       |
|                        | IDs        |                    |included in the network      | 
|                        |            |                    |construction.                |
+------------------------+------------+--------------------+-----------------------------+ 
|Threshold               | Integer    |``--threshold``     |Distance threshold when      | 
|                        |            |                    |constructing network.        |
+------------------------+------------+--------------------+-----------------------------+ 
|Step                    | Integer    |``--step``          |Step to use when iterating   | 
|                        |            |                    |through trajectory frames.   |
+------------------------+------------+--------------------+-----------------------------+ 
|Generate plots          | Boolean    |``--generate-plots``|Set to generate figures.     | 
+------------------------+------------+--------------------+-----------------------------+ 
|Calculate BC            | Boolean    |``--calc-BC``       |Set to calculate average     | 
|                        |            |                    |shortest path matrix for the | 
|                        |            |                    |network                      |
+------------------------+------------+--------------------+-----------------------------+
|Calculate L             | Boolean    |``--calc-L``        |Set to calculate betweenness | 
|                        |            |                    |centrality matrix for the    | 
|                        |            |                    |network                      |
+------------------------+------------+--------------------+-----------------------------+ 
|Discard graphs          | Boolean    |``--discard-graphs``|Set to discard the network   | 
|                        |            |                    |once BC and L matrices have  | 
|                        |            |                    |been calculated              |
+------------------------+------------+--------------------+-----------------------------+ 

Given a trajectory called ``minimized.dcd`` and a topology file called ``minimized.pdb``, the following command could be used: ::

	python calc_network.py --topology minimized.pdb --threshold 7.0 --step 100 --generate-plots --calc-BC --calc-L --discard-graphs minimized.dcd

The above command will calculate the network for every 100th frame in the trajectory. Depending on the size of your trajectory, you may want to increase this ``--step``. Edges in the network will be created between nodes that are within 7 Angstroms of each other. The average shortest path for each residue in each frame and the betweenness centrality of each residue in each frame will be calculated as both flags have been set in the above command. In addition, the ``--discard-graphs`` flag was set. As such, the networks for each frame will be discarded once BC and L have been calculated, saving disk space. By default, the networks for each frame are save in both ``gml`` and ``graphml`` format.

**Outputs:**

================  ===================================================================================================================================================================
Output            Description
================  ===================================================================================================================================================================
BC Matrices       For each frame analyzed, an Nx1 matrix is produced, where N is the number of residues in the protein and each value represents the BC for the residue at that index
avg_L Matrices    For each frame analyzed, an Nx1 matrix is produced, where N is the number of residues in the protein and each value represents the L to the residue at that index
BC & L Plots      If ``--generate-plots`` flag is set, PNG figures are produced for the BC and L matrices
Network graphs    If ``--discard-graphs`` flag is set, do not save the networks produced for each frame
================  ===================================================================================================================================================================

Calculating ΔL
----------------------

If the ``--calc-L`` flag in the previous command is set, a number of Nx1 L matrices will be generated. Given the trajectory ``minimized.dcd``, the matrices will be named ``minimized_<frame>_avg_L.dat``, where ``<frame>`` is the frame index in the trajectory. 

**Command:** :: 
	
	python calc_delta_L.py <options> --reference <frame> --alternatives <frames>

**Inputs:**

=========================  ===========  ====================  ========================================================================================================================================================
 Input (*\*required*)      Input type   Flag                  Description                  
=========================  ===========  ====================  ========================================================================================================================================================
Reference frame *          File         ``--reference--``     Nx1 matrix to be used as the reference (normally the frame from time 0). Delta L will be worked out by comparing the alternative frames to this one.    
Alternative frames *       File/s       ``--alternatives``    The remaining Nx1 matrices that should be compared to the reference matrix
Normalize                  Boolean      ``--normalize``       Set this flag to normalize the values (ΔL/L)
Generate plots             Boolean      ``--generate-plots``  Set to generate figures
=========================  ===========  ====================  ========================================================================================================================================================

Given a set of average shortest path .dat files ``minimized_*_avg_L.dat`` (generated with ``calc_network.py``), the ``minimized_0_avg_L.dat`` file could be used as the reference and the rest could be used as the alternatives. If the reference .dat file is renamed to ``reference.dat``, the following command could be used: ::

	python calc_delta_L.py --normalize --generate-plots --reference reference.dat --alternatives minimized_*_avg_L.dat

The above command will generate plots as well as Nx1 matrices representing the difference in L between each alternative and the reference frame. The values will be normalized by dividing by the reference values (ΔL/L).

**Outputs:**

================  ===================================================================================================================================================================
Output            Description
================  ===================================================================================================================================================================
ΔL Matrices       Nx1 matrices representing the change in L between the reference matrix and each alternative
ΔL Plots          Figures for each alternative frame, plotting the difference between L in the alternative and reference
================  ===================================================================================================================================================================

Calculating ΔBC
-----------------------

If the ``--calc-BC`` flag was set when running the ``calc_network.py`` script, a number of Nx1 BC matrices will be generated. Given the trajectory ``minimized.dcd``, the matrices will be named ``minimized_<frame>_bc.dat``, where ``<frame>`` is the frame index in the trajectory. 

**Command:** :: 
	
	python calc_delta_BC.py <options> --reference <frame> --alternatives <frames>

**Inputs:**

=========================  ===========  ====================  ========================================================================================================================================================
 Input (*\*required*)      Input type   Flag                  Description                  
=========================  ===========  ====================  ========================================================================================================================================================
Reference frame *          File         ``--reference--``     Nx1 matrix to be used as the reference (normally the frame from time 0). Delta BC will be worked out by comparing the alternative frames to this one.    
Alternative frames *       File/s       ``--alternatives``    The remaining Nx1 matrices that should be compared to the reference matrix
Generate plots             Boolean      ``--generate-plots``  Set to generate figures
=========================  ===========  ====================  ========================================================================================================================================================

Given a set of BC .dat files ``minimized_*_bc.dat`` (generated with ``calc_network.py``), the ``minimized_0_bc.dat`` file could be used as the reference and the rest could be used as the alternatives. If the reference .dat file is renamed to ``reference.dat``, the following command could be used: ::

	python calc_delta_BC.py --generate-plots --reference reference.dat --alternatives minimized_*_bc.dat

The above command will generate plots as well as Nx1 matrices representing the difference in BC between each alternative and the reference frame.

**Outputs:**

================  ===================================================================================================================================================================
Output            Description
================  ===================================================================================================================================================================
ΔBC Matrices      Nx1 matrices representing the change in BC between the reference matrix and each alternative
ΔBC Plots         Figures for each alternative frame, plotting the difference between BC in the alternative and reference
================  ===================================================================================================================================================================


Calculating Average BC and L (and standard deviation)
-----------------------------------------------------

The ``avg_network.py`` script can be used to calculate and plot the average BC and L as well as the standard deviation of these measurements over the course of the trajectory.

**Command:** ::
	
	python avg_network.py <options> --data-type <BC/delta-BC/L/delta-L> --data <matrices>

**Inputs:**

=========================  ===========  ====================  ========================================================================================================================================================
 Input (*\*required*)      Input type   Flag                  Description                  
=========================  ===========  ====================  ========================================================================================================================================================    
Data *                     File/s       ``--data``            The .dat files that will be averaged 
Data types *               Text         ``--data-type``       Type of data - BC/delta-BC/L/delta-L
Prefix                     Text         ``--prefix``          Prefix used to name outputs
Generate plots             Boolean      ``--generate-plots``  Generate figures/plots     
X axis label               Text         ``--x-label``         Label for x-axis (use $\Delta$ for delta sign)
Y axis label               Text         ``--y-label``         Label for y-axis (use $\Delta$ for delta sign)
Max Y axis value           Integer      ``--y-max``           Maximum value on y-axis
Min Y axis value           Integer      ``--y-min``           Minimum value on y-axis
Graph title                Text         ``--title``           Title of plot (use $\Delta$ for delta sign)
X-axis start value         Integer      ``--initial-x``       The start index of the X-axis
Split position             Integer      ``--split-pos``       Position to split the network at for large networks. Splits the plot at the given position to create two plots. Useful when analysing a dimer.    
Graph title 1              Text         ``--title-1``         Title of first plot  
Graph title 2              Text         ``--title-2``         Title of second plot  
X-axis start value 1       Integer      ``--initial-x-1``     The start index of the x-axis for the first plot      
X-axis start value 2       Integer      ``--initial-x-2``     The start index of the x-axis for the second plot                
=========================  ===========  ====================  ========================================================================================================================================================

Given a set of .dat files generated by one of the previous commands (e.g. ``minimized_*_bc_delta_BC.dat``), the following command could be used: ::
	
	python avg_network.py --data minimized_*_bc_delta_BC.dat --data-type delta-BC --prefix minimized --generate-plots --x-label "Residues" --y-label "Avg delta BC" --title "My Protein"

The above command will generate two new .dat files and a PNG plot. The first .dat file, ``minimized_avg.dat``, contains an Nx1 matrix with the average ΔBC values for each residue over the course of the simulation. The second .dat file, ``minimized_std_dev.dat``, contains the standard deviation of ΔBC for each residue over the course of the simulation. The graph plots residues on the X axis and ΔBC on the Y axis. The average values are shown as a line and the standard deviation, representing the fluctuation of ΔBC over the course of the trajectory, are shown as error bars over each residue.

**Outputs:**

=================  ===================================================================================================================================================================
Output             Description
=================  ===================================================================================================================================================================
Average .dat file  Nx1 matrix representing the average BC/ΔBC/L/ΔL values from the input matrics
Std dev .dat file  Nx1 matrix representing the standard deviation of the BC/ΔBC/L/ΔL values of the input matrics 
Plot               The plotted values from the above matrices 
=================  =================================================================================================================================================================== 

SNP Analysis - comparing mutant to wild type trajectories
---------------------------------------------------------



SNP Analysis - residue contact maps
---------------------------------------------------------





