Installation
========================================

Platform compatibility
-------------------------------

MD-TASK should be compatible with any Linux/Unix-based platform, although installation of system dependencies may differ. It has been successfully tested on the following platforms:

- Ubuntu Linux
- MacOS
- Windows 10 (with bash)

Install system dependencies
-----------------------------

*Note: package version numbers may differ depending on the OS version. For example, in Ubuntu 16.04, `libpng12-dev` must be installed. However, in Ubuntu 17.04, `libpng-dev` should be installed.*

**Ubuntu 16.04:** ::

	sudo apt-get install virtualenvwrapper python-dev libblas-dev liblapack-dev libatlas-base-dev gfortran libpng12-dev libfreetype6-dev python-tk r-base

**Windows 10:** 

1. Enable the Windows Subsystem for Linux (WSL) by following `these instructions <https://msdn.microsoft.com/en-us/commandline/wsl/install_guide>`_.

2. Install the system dependencies as with Ubuntu above.

**MacOS:**

1. On MacOS, Python comes installed by default, but the default version my not be ideal. Follow `these instructions <http://exponential.io/blog/2015/02/11/install-python-on-mac-os-x-for-development/>`_ to install a more up-to-date version of Python.

2. Next, install virtualenv by following `these instructions <http://exponential.io/blog/2015/02/10/install-virtualenv-and-virtualenvwrapper-on-mac-os-x/>`_

Install Python dependencies
--------------------------------

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

Download the project
-------------------------------

MD-TASK can be cloned from it's GitHub repository ::

	git clone https://github.com/RUBi-ZA/MD-TASK.git
	cd MD-TASK

Always activate the virtual environment you created in the previous step when using MD-TASK.
