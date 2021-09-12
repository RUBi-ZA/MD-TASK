#!/bin/bash

python -m unittest discover
./test_CM.sh
./test_DCC.sh
./test_PRS.sh
./test_CP.sh
