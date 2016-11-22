Usage
=========


Network Analysis
-----------------

**Measurements**

*1. Betweenness Centrality (BC)*

Betweenness centrality (BC) is a measure of how important a residue is for communication within a protein. It is equal to the number of shortest paths from all vertices to all others that pass through that node. Residues in a protein that have a high BC reveal locations that tend to be important for controlling inter-domain commuication in a protein.

*2. Average Shortest Path(L)*

The average shortest path (L) to a given residue is calculated by working out all the shortest paths to the given node and dividing by the number of paths. The average shortest path to a residue gives an idea of how accessible the residue is within the protein. This can be used to, for example, analyze SNPs. A mutation may result in a change in L of a number of residues in the protein. This may indicate that the mutation has an important effect on protein function e.g. previous studies have suggested that positions that result in high delta L values may steer conformational changes.

**Calculating the network**

Required inputs:

+------------------------+------------+----------+----------+
| Input                  | Input type | Description         |
| (header rows optional) |            |                     |
+========================+============+=====================+
|                        |            |                     |
+------------------------+------------+---------------------+
|                        |            |                     |
+------------------------+------------+---------------------+  

=====  =====  =======
A      B      A and B
=====  =====  =======
False  False  False
True   False  False
False  True   False
True   True   True
=====  =====  =======