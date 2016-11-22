# MD-RIN
Tool suite for analysing molecular dynamics trajectories using network analysis and PRS

##1. Introduction

MD-RIN consists of a suite of Python scripts that have been developed to analyze molecular dynamics trajectories. These scripts fall into 3 categories:

###1.1. Network analysis

Residue Interaction Networks (RIN) are analyzed using a branch of Mathematics known as graph theory. In a RIN, each residue in the protein is a node in the network. An edge (or connection) between two nodes exists if there is an interaction between the two residues that they represent. In RINs, an interaction between two residues exists if the residues are within a user-defined cut-off distance (usually around 6.5 – 7.5 Å) of each other. Once the network has been constructed, there are various network measures that can be used to analyze it. Currently, MD-RIN can be used to analyze the change in betweenness centrality (BC) and average shortest path (L) of residues in a protein over a molecular dynamics simulation.

###1.1.1. Betweenness Centrality

Betweenness centrality (BC) is a measure of how important a residue is for communication within a protein. It is equal to the number of shortest paths from all vertices to all others that pass through that node. Residues in a protein that have a high BC reveal locations that tend to be important for controlling inter-domain commuication in a protein. 

###1.1.2. Average Shortest Path

The average shortest path (L) to a given residue is calculated by working out all the shortest paths to the given node and dividing by the number of paths. The average shortest path to a residue gives an idea of how accessible the residue is within the protein. This can be used to, for example, analyze SNPs. A mutation may result in a change in L of a number of residues in the protein. This may indicate that the mutation has an important effect on protein function e.g. previous studies have suggested that positions that result in high delta L values may steer conformational changes.

###2.2. Perturbation Response Scanning

###2.3. MD Correlation

##2. Installation

###2.1. Download the project
```bash
git clone https://github.com/RUBi-ZA/JMS.git
cd MD-RIN
```
###2.2. Install dependencies and set up Python virtual environment
```bash
sudo apt-get install virtualenvwrapper python-dev libblas-dev liblapack-dev libatlas-base-dev gfortran libpng12-dev libfreetype6-dev python-tk
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy 
pip install scipy 
pip install matplotlib cython networkx natsort
pip install mdtraj
```


