# MDM-TASK-web [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

Web server tools for MD-TASK and MODE-TASK.

Take a look at our web application [MDM-TASK-web](https://mdmtaskweb.rubi.ru.ac.za/) and find the [BioRxiv article here](https://www.biorxiv.org/content/10.1101/2021.01.29.428734v1)


## Installation
1. Clone the MDM-TASK-web branch
```
git clone -b mdm-task-web https://github.com/RUBi-ZA/MD-TASK
```
Go inside the cloned MD-TASK folder and run the following to create the working environment named "mdmtaskweb"
```
conda env create -f environment.yml
```
Activate the conda environment
```
conda activate mdmtaskweb
```
The scripts are located in the root directory of the project. To run the scripts, ensure that the conda environment is activated. 

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RUBi-ZA/MD-TASK/mdmtask-dev?filepath=example%2Fmdmtaskweb_tutorial.ipynb)
