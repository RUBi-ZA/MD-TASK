#!/usr/bin/env python
"""
Different engines for coarse-graining trajectories
"""
import re
import os
import warnings
from shlex import quote as q
from subprocess import STDOUT, DEVNULL
from subprocess import check_output as call
from tempfile import NamedTemporaryFile
from lib.libcheck import _check_filename

def _not_implemented(method=None):
    """Checks implementation"""
    if method is not None:
        warnings.warn("{} not implemented".format(method))
        return NotImplemented
    return None

def basename(text, complete=False):
    """Strips the file extension.
    complete strips the directory"""
    if complete:
        text = os.path.basename(text)
    return os.path.splitext(text)[0]

def reduce_mdanalysis(topology_fname, trajectory_fname, selection='name CA CB'):
    """Save only CA and CB atoms"""
    import MDAnalysis as mda
    topology_out = "{}_small.pdb".format(basename(topology_fname, complete=True))
    trajectory_out = "{}_small.xtc".format(basename(trajectory_fname,
                                                    complete=True))
    reduced_trajectory = mda.Universe(topology_fname,
            trajectory_fname).select_atoms(selection)
    reduced_topology = mda.Universe(topology_fname).select_atoms(selection)
    reduced_topology.write(topology_out)
    print("INFO: Wrote {}".format(topology_out))
    reduced_trajectory.write(trajectory_out, frames="all")
    print("INFO: Wrote {}".format(trajectory_out))

def reduce_mdtraj(topology_fname, trajectory_fname, selection="name CA CB"):
    """Save only CA and CB atoms"""
    import mdtraj as md
    topology_out = "{}_small.pdb".format(basename(topology_fname, complete=True))
    trajectory_out = "{}_small.xtc".format(basename(trajectory_fname, complete=True))
    topology = md.load(topology_fname)
    topology.atom_slice((topology.top.select(selection)), inplace=True)
    trajectory = md.load(trajectory_fname, top=topology_fname)
    trajectory.atom_slice(trajectory.top.select(selection), inplace=True)
    topology.save_pdb(topology_out)
    print("INFO: Wrote {}".format(topology_out))
    trajectory.save_xtc(trajectory_out)
    print("INFO: Wrote {}".format(trajectory_out))

def reduce_pytraj(topology_fname, trajectory_fname, selection="@CA,CB"):
    """Save only CA and CB atoms"""
    import pytraj as pt
    topology_out = "{}_small.pdb".format(basename(topology_fname, complete=True))
    trajectory_out = "{}_small.xtc".format(basename(trajectory_fname,
                                                    complete=True))
    trajectory = pt.load(trajectory_fname, top=topology_fname)
    topology = pt.load(topology_fname)
    topology = topology[selection]
    trajectory = trajectory[selection]
    topology.save(topology_out)
    print("INFO: Wrote {}".format(topology_out))
    trajectory.save(trajectory_out)
    print("INFO: Wrote {}".format(trajectory_out))

def reduce_gmx(topology_fname, trajectory_fname, selection=None):
    """Save only CA and CB atoms"""
    _not_implemented(selection)
    topology_out = "{}_small.pdb".format(basename(topology_fname, complete=True))
    trajectory_out = "{}_small.xtc".format(basename(trajectory_fname,
                                                    complete=True))
    index_out = "{}_small.ndx".format(basename(topology_fname, complete=True))
    topology_base = basename(topology_fname)
    trajectory_base = basename(trajectory_fname)
    _check_filename(topology_fname)
    _check_filename(trajectory_fname)
    makeindex_cmd = "echo 'a CA CB\\nq' | gmx make_ndx -f {topology_fname} -o {index_out}".format(
        topology_fname=q(topology_fname), index_out=q(index_out))
    topology_cmd = "echo 'CA_CB' | gmx trjconv -s {topology_fname} -f {topology_fname} -n {index_out} -o {topology_out}".format(
        topology_fname=q(topology_fname), index_out=q(index_out),
        topology_out=q(topology_out))
    trajectory_cmd = "echo 'CA_CB' | gmx trjconv -s {topology_fname} -f {trajectory_fname} -n {index_out} -o {trajectory_out}".format(
        topology_fname=q(topology_fname), trajectory_fname=q(trajectory_fname),
        index_out=q(index_out), trajectory_out=q(trajectory_out))
    call(makeindex_cmd, shell=True, stderr=DEVNULL)
    call(topology_cmd, shell=True, stderr=DEVNULL)
    print("INFO: Wrote {}_small.pdb".format(topology_base))
    call(trajectory_cmd, shell=True, stderr=DEVNULL)
    print("INFO: Wrote {}_small.xtc".format(trajectory_base))

def reduce_cpptraj(topology_fname, trajectory_fname,
               selection="!@CA,CB"):
    """Save only CA and CB atoms"""
    topology_out = "{}_small.pdb".format(basename(topology_fname, complete=True))
    trajectory_out = "{}_small.xtc".format(basename(trajectory_fname, complete=True))
    command1 = "echo 'parm {topology_fname}\n\
                trajin {topology_fname}\n\
                strip {selection}\n\
                trajout {topology_out}\ngo\nquit' | cpptraj".format(
                        topology_fname=q(topology_fname),
                        selection=q(selection),
                        topology_out=q(topology_out))
    command2 = "echo 'parm {topology_fname}\n\
                trajin {trajectory_fname}\n\
                strip {selection}\n\
                trajout {trajectory_out}\ngo\nquit' | cpptraj".format(
                        topology_fname=q(topology_fname),
                        trajectory_fname=q(trajectory_fname),
                        selection=q(selection),
                        trajectory_out=q(trajectory_out))
    call(command1, shell=True)
    print("INFO: Wrote {}".format(topology_out))
    call(command2, shell=True)
    print("INFO: Wrote {}".format(trajectory_out))

def reduce_vmd(topology_fname, trajectory_fname,
               selection="(name CA or name CB) and not solvent"):
    """Save only CA and CB atoms, adapted from
    http://md-task.readthedocs.io/en/latest/general.html
    """
    topology_out = "{}_small.pdb".format(
            basename(topology_fname, complete=True))
    trajectory_out = "{}_small.dcd".format(
            basename(trajectory_fname, complete=True))
    trajext = os.path.splitext(trajectory_fname)[1][1:]
    #Selection
    commands = 'mol new {topfname}\n\
                set s1 [atomselect top "{selection}"]\n\
                animate write pdb {topout} sel $s1\n\
                animate read {trajext} {trajfname} waitfor all\n\
                animate write dcd {trajout} waitfor all sel $s1\n\
                quit'.format(topfname=q(topology_fname),
                             selection=q(selection),
                             topout=q(topology_out),
                             trajext=q(trajext),
                             trajfname=q(trajectory_fname),
                             trajout=q(trajectory_out))
    tmpfile = NamedTemporaryFile('wt')
    tmpfile.file.write(commands)
    tmpfile.file.close()
    output = call('vmd -dispdev none -e {}'.format(q(tmpfile.name)),
                  shell=True, stderr=STDOUT)
    vmderror = "\n".join(re.findall('ERROR.+', output.decode('utf8')))
    if vmderror != "":
        print(vmderror)
        try:
            os.remove(topology_out)
            os.remove(trajectory_out)
        except IOError:
            pass
    else:
        print("INFO: Wrote {}".format(topology_out))
        print("INFO: Wrote {}".format(trajectory_out))


engines = {
        "MDAnalysis": reduce_mdanalysis,
        "pytraj": reduce_pytraj,
        "mdtraj": reduce_mdtraj,
        "gmx": reduce_gmx,
        "cpptraj": reduce_cpptraj,
        "vmd": reduce_vmd
        }
