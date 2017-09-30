#!/bin/bash

virtualenv venv
source venv/bin/activate

pip install --upgrade pip
pip install numpy==1.13.3
pip install scipy==0.19.1
pip install matplotlib==2.0.2
pip install cython==0.27 
pip install networkx==1.11
pip install natsort==5.1.0
pip install pandas==0.20.3
pip install mdtraj==1.9.1
