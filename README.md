# MD-RIN
Tool suite for analysing molecular dynamics trajectories using network analysis and PRS

## Installation

###1. Download the project

```bash
git clone https://github.com/RUBi-ZA/JMS.git
cd MD-RIN
```

###2. Install dependencies and set up Python virtual environment

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

## Usage

MD-RIN consists of a suite of Python scripts that have been designed to analyze molecular dynamics simulations. These scripts fall into 3 categories:

###1. Network analysis

###2. Perturbation Response Scanning

###3. MD Correlation
