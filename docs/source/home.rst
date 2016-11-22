MD-RIN
==========

MD-RIN consists of a suite of Python scripts that have been developed to analyze molecular dynamics trajectories. These scripts fall into 3 categories:

**1. Residue Interaction Network (RIN) analysis**

Residue Interaction Networks (RIN) are analyzed using a branch of Mathematics known as graph theory. In a RIN, each residue in the protein is a node in the network. An edge (or connection) between two nodes exists if there is an interaction between the two residues that they represent. In RINs, an interaction between two residues exists if the residues are within a user-defined cut-off distance (usually around 6.5 – 7.5 Å) of each other. Once the network has been constructed, there are various network measures that can be used to analyze it. Currently, MD-RIN can be used to analyze the change in betweenness centrality (BC) and average shortest path (L) of residues in a protein over a molecular dynamics simulation.


**2. Perturbation Response Scanning (PRS)**


**3. MD Correlation calculation**



Contribute
----------------

- Issue Tracker: https://github.com/RUBi-ZA/MD-RIN/issues
- Source Code: https://github.com/RUBi-ZA/MD-RIN


License
---------------

The project is licensed under GNU GPL 3.0