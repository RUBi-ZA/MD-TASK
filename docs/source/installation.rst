Installation
========================================

Platform compatibility
-------------------------------

MD-TASK should be compatible with any Linux/Unix-based platform, although installation of system dependencies may differ. It has been successfully tested on the following platforms:

- Ubuntu Linux
- MacOS
- Windows 10 (with bash)

Download the project
-------------------------------

MD-TASK can be cloned from it's GitHub repository ::

	git clone https://github.com/RUBi-ZA/MD-TASK.git
	cd MD-TASK

Install dependencies
---------------------

Install system dependencies (e.g. Ubuntu) ::

	sudo apt-get install virtualenvwrapper python-dev libblas-dev liblapack-dev libatlas-base-dev gfortran libpng12-dev libfreetype6-dev python-tk r-base


We recommend using a Python virtual environment when using MD-TASK ::

	virtualenv venv
	source venv/bin/activate
	pip install --upgrade pip
	pip install numpy 
	pip install scipy 
	pip install matplotlib cython networkx natsort
	pip install mdtraj


Install the igraph package for R: ::

	R
	> install.packages("igraph")
