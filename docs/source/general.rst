General
========

Add MD-RIN to your PATH
-------------------------

To make tools in the MD-RIN suite available from anywhere on the command line, add the root MD-RIN directory to your PATH environment variable as follows: ::

	export PATH=/path/to/MD-RIN:$PATH


Trajectory vs Topology
------------------------

Most of the MD-RIN tools require both a trajectory file and topology/structure file as input. This is because most trajectory formats only contain the atom co-ordinates and not the topological information such as atom and residue names, chains, and bond information. The topology file can be the PDB file that was used in the molecular dynamics simulation to produce the trajectory. When supplying these files, it is important to note that the trajectory file and topology file must contain the exact same number of atoms i.e. if the trajectory has been reduced to CA and CB atoms only (as described below), the topology file must be reduced to the same.

Reducing your trajectory
-----------------------------

Molecular dynamics trajectories can be extremely large. However, MD-RIN tools only require the alpha and beta carbon atoms to be present. To save space and improve performance, the following VMD script can be used to reduce a trajectory, to the bare essentials: ::

	mol load pdb example.pdb
	set s1 [atomselect top "name CA or name CB and not solvent"]
	animate write pdb example.pdb sel $s1 
	animate read xtc example.xtc waitfor all
	animate write dcd example_small.dcd waitfor all sel $s1 
	quit

The above assumes that your topology file is a PDB file named ``example.pdb`` and that your trajectory is named ``example.xtc``. It then writes out the reduced structure and trajectory to ``example.pdb` and ``example_small.dcd`` respectively. You should change these names accordingly. You will also note that the above converts the trajectory from XTC to DCD format. This is not necessary. 

For very large trajectories that do not fit in memory, reducing as shown above may be necessary. Note that when reducing the trajectory, the same reduction should be applied to the topology PDB file i.e. the trajectory and topology files should have the exact same number of atoms.

Test Data
----------

There is test data located in the 'examples' directory. Four files are included here:

=====================  =======================================================================================================================================================
File                    Description
=====================  =======================================================================================================================================================
``example_small.dcd``  An example trajectory that has been reduced to alpha and beta carbons only
``example.pdb``        A PDB file that corresponds to the above trajectory (to be used for topology information)
``Initial.txt``        An XYZ co-ordinate file representing the initial conformation of a protein (used for PRS)
``Final.txt``          An XYZ co-ordinate file representing the target conformation of a protein (used for PRS)
=====================  =======================================================================================================================================================

Logging
--------

All scripts in the suite have two arguments for logging. By default, logging is switched on and is written to the terminal. This can be changed with the following arguments:

============  ==================  =====================================================================================================================================
Input         Flag                Description
============  ==================  =====================================================================================================================================
Log file      ``--log-file``      Provide a path to a file that will store the logging output of the command. By default, the output will be written to the terminal.
Silent        ``--silent``        Switch off logging
============  ==================  =====================================================================================================================================
