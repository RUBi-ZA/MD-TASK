# MD-RIN
Tool suite for analysing molecular dynamics trajectories using network analysis and PRS

##1. Installation

###1.1. Download the project

```bash
git clone https://github.com/RUBi-ZA/JMS.git
cd MD-RIN
```

###1.2. Install dependencies and set up Python virtual environment

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

##2. Usage

MD-RIN consists of a suite of Python scripts that have been designed to analyze molecular dynamics simulations. These scripts fall into 3 categories:

###2.1. Network analysis

Residue Interaction Networks (RIN) are analyzed using a branch of Mathematics known as graph theory. In a RIN, each residue in the protein is a node in the network. An edge (or connection) between two nodes exists if there is an interaction between the two residues that they represent. In RINs, an interaction between two residues exists if the residues are within a user-defined cut-off distance (usually around 6.5 – 7.5 Å) of each other. Once the network has been constructed, there are various network measures that can be used to analyze it. Currently, MD-RIN can be used to analyze the change in betweenness centrality (BC) and average shortest path (L) of residues in a protein over a molecular dynamics simulation.

###2.1.1. Betweenness Centrality

Betweenness centrality (BC) is a measure of how important a residue is for communication within a protein. It is equal to the number of shortest paths from all vertices to all others that pass through that node.

###2.1.2. Average Shortest Path

The average shortest path (L) to a given residue is calculated by working out all the shortest paths to the given node and dividing by the number of paths. The average shortest path to a residue gives and idea of how accessible the residue is within the protein.

###2.2. Perturbation Response Scanning

###2.3. MD Correlation
