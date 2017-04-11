#!/bin/bash

virtualenv venv
source venv/bin/activate

pip install --upgrade pip
pip install numpy
pip install scipy
pip install matplotlib cython networkx natsort
pip install mdtraj
