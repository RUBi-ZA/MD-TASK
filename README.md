[![Build Status](https://travis-ci.org/RUBi-ZA/MD-TASK.svg?branch=master)](https://travis-ci.org/RUBi-ZA/MD-TASK)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/ece0a51e7cf4436abf71795051f4ee7b)](https://www.codacy.com/gh/RUBi-ZA/MD-TASK?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=RUBi-ZA/MD-TASK&amp;utm_campaign=Badge_Grade)

# MD-TASK

A suite of Python scripts that have been developed to analyze molecular dynamics trajectories. Detailed documentation on how to install and use MD-TASK can be found on our [ReadTheDocs](http://md-task.readthedocs.io/en/latest/index.html) site.

## Installation

*1. Download the project:*
```bash
git clone  https://github.com/RUBi-ZA/MD-TASK.git
cd MD-TASK
```

*2. Create and activate the `mdtask` environment*

You can either run:

```bash
conda create -n mdtask pandas conda-forge::matplotlib conda-forge::natsort anaconda::networkx conda-forge::mdtraj
conda activate mdtask
```

or you can use the provided YAML file:

```bash
conda env create -f md-task.yml
conda activate md-task
```

## Usage
The scripts are located in the root directory of the project. To run the scripts, ensure that the conda virtual environment is activated. For more information, please refer to the detailed documentation [here](http://md-task.readthedocs.io/en/latest/index.html).
