# MDM-TASK-web [![fair-software.eu](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B%20%20%E2%97%8F%20%20%E2%97%8B-orange)](https://fair-software.eu)

**MDM-TASK-web** is web server for MD-TASK and MODE-TASK and can be accessed here: [https://mdmtaskweb.rubi.ru.ac.za](https://mdmtaskweb.rubi.ru.ac.za)

If you use MDM-TASK-web in your work, please cite the following:
 - [MDM-TASK-web: MD-TASK and MODE-TASK web server for analyzing protein dynamics](https://doi.org/10.1016/j.csbj.2021.08.043)
 - [MD-TASK: a software suite for analyzing molecular dynamics trajectories](https://dx.doi.org/10.1093%2Fbioinformatics%2Fbtx349)

## Installation
Clone the MDM-TASK-web branch
```
git clone -b mdm-task-web  https://github.com/RUBi-ZA/MD-TASK
```
Go inside the cloned MD-TASK folder and run the following to create the working environment named "mdmtaskweb"
```
conda env create -f environment.yml
```
Activate the conda environment
```
conda activate mdmtaskweb
```
 - Scripts are located in the src directory of the project. To run the scripts, ensure that the conda environment is activated. 
 - Example topology and (short) trajectory files are in the data folder
 - Click on the binder logo for a tutorial of the command line
 
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/RUBi-ZA/MD-TASK/mdmtask-dev?filepath=tutorial%2Fmdmtaskweb_tutorial.ipynb)
