<img src="https://travis-ci.org/RUBi-ZA/MD-TASK.svg?branch=master" align="right">

# MD-TASK

Tool suite for analysing molecular dynamics trajectories using network analysis and PRS.

## Introduction

MD-TASK consists of a suite of Python scripts that have been developed to analyze molecular dynamics trajectories. Detailed documentation can be found on our [ReadTheDocs](http://md-task.readthedocs.io/en/latest/index.html) site.

## Installation

*Download the project:*
```bash
git clone https://github.com/RUBi-ZA/MD-TASK.git
cd MD-TASK
```
*Install dependencies and set up Python virtual environment:*
```bash
sudo apt-get install virtualenvwrapper python-dev libblas-dev liblapack-dev libatlas-base-dev gfortran libpng12-dev libfreetype6-dev python-tk r-base
virtualenv venv
source venv/bin/activate
pip install --upgrade pip
pip install numpy
pip install scipy
pip install matplotlib cython networkx natsort
pip install mdtraj
```
*Install igraph package for R:*
```bash
R
> install.packages("igraph")
```

## Usage

The scripts are located in the root directory of the project. To run the scripts, ensure that the virtual environment is activated. For more info, find detailed documentation [here](http://md-task.readthedocs.io/en/latest/index.html). 

