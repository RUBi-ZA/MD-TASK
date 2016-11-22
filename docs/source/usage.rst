Network Analysis
==================

Residue Interaction Networks (RIN) are analyzed using a branch of Mathematics known as graph theory. In a RIN, each residue in the protein is a node in the network. An edge (or connection) between two nodes exists if there is an interaction between the two residues that they represent. In RINs, an interaction between two residues exists if the residues are within a user-defined cut-off distance (usually around 6.5 – 7.5 Å) of each other. Once the network has been constructed, there are various network measures that can be used to analyze it. Currently, MD-RIN can be used to analyze the change in betweenness centrality (BC) and average shortest path (L) of residues in a protein over a molecular dynamics simulation.

Measurements
-----------------

*1. Betweenness Centrality (BC)*

Betweenness centrality (BC) is a measure of how important a residue is for communication within a protein. It is equal to the number of shortest paths from all vertices to all others that pass through that node. Residues in a protein that have a high BC reveal locations that tend to be important for controlling inter-domain commuication in a protein.

*2. Average Shortest Path(L)*

The average shortest path (L) to a given residue is calculated by working out all the shortest paths to the given node and dividing by the number of paths. The average shortest path to a residue gives an idea of how accessible the residue is within the protein. This can be used to, for example, analyze SNPs. A mutation may result in a change in L of a number of residues in the protein. This may indicate that the mutation has an important effect on protein function e.g. previous studies have suggested that positions that result in high delta L values may steer conformational changes.

Calculating the network
------------------------

*Command:* ::
	
	python calc_network <options> <trajectory>

*Inputs:*

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

*Outputs:*

For each frame in the the trajectory that a RIN was generated for, an 2 Nx1 matrices will be printed to file, one for BC and one for L (if both flags were set), where N is the number of residues in the protein. Each value in these matrices represents BC or L for the residue at that index.

Two PNG plots will also be generated for each frame if the ``--generate-plots`` flag was set.
