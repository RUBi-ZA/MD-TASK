#! /bin/bash

python3 -m venv venv
. venv/bin/activate

pip3 install --upgrade pip
pip3 install numpy==1.18.1
pip3 install scipy==1.4.1
pip3 install matplotlib==3.2.0
pip3 install cython==0.29.15 
pip3 install networkx==1.11
pip3 install natsort==7.0.1
pip3 install pandas==0.25.3
pip3 install mdtraj==1.9.3
pip3 install seaborn==0.9.1
