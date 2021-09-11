"""Modules for general checks"""
import os
import sys
import importlib
from subprocess import call, DEVNULL

modules = ["MDAnalysis", "mdtraj", "pytraj"]
tools = ["gmx", "cpptraj", "vmd"]

def _check_filename(filename):
    """Checks if argument is a path or file"""
    if os.path.isdir(filename):
        print("ERROR: Path {} is a directory".format(filename))
        sys.exit(1)
    if not os.path.exists(filename):
        print("ERROR: File {} does not exist".format(filename))
        sys.exit(1)

def _check_module():
    """Checks which module is available"""
    for module in modules:
        try:
            globals()[module] = importlib.import_module(module)
            print("INFO: Using the {0} module.".format(module))
            return module
        except ModuleNotFoundError:
            print("WARN: Module {0} was not found.".format(module.split()[0]))
    return False

def _check_tool():
    """Check which tool is available"""
    for tool in tools:
        if tool == "vmd":
            state = call("echo|{} -dispdev none".format(tool),
                         shell=True, stderr=DEVNULL,
                         stdout=DEVNULL)
        else:
            state = call("echo|{}".format(tool),
                         shell=True, stderr=DEVNULL,
                         stdout=DEVNULL)
        if state == 0:
            print("INFO: Using {0}.".format(tool))
            return tool
        else:
            print("WARN:{0} was not found.".format(module))
    return False

def check_tools():
    """Checks and sets the library/tool to use"""
    print("INFO: Autodetecting engine.")
    module = _check_module()
    if module:
        return module
    else:
        tool = _check_tool()
        if tool:
            return tool
