import sys

def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)

class LogLevel(object):
    DEBUG = 1
    INFO = 2
    WARN = 3
    ERROR = 4

class Logger(object):

    def __init__(self, stream=sys.stdout, silent=False, log_level=LogLevel.INFO):
        self.stream = stream
        self.silent = silent
        self.log_level = log_level

    def _log(self, message):
        if not self.silent:
            self.stream.write(message)
            self.stream.flush()

    def debug(self, message):
        if self.log_level <= LogLevel.DEBUG:
            message = "DEBUG::%s" % message
            self._log(message)

    def info(self, message):
        if self.log_level <= LogLevel.INFO:
            message = "INFO::%s" % message
            self._log(message)

    def warn(self, message):
        if self.log_level <= LogLevel.WARN:
            message = "WARNING::%s" % message
            self._log(message)

    def error(self, message):
        if self.log_level <= LogLevel.ERROR:
            message = "ERROR::%s" % message
            self._log(message)
def dat2xmgrace(val, prefix, output, traj, selection="(name CB and protein) or (name CA and resname GLY)"):
    import mdtraj as md
    import pandas as pd
    import os

    assert type(traj) is md.core.trajectory.Trajectory, "traj has to be an object of type md.core.trajectory.Trajectory"

    protein = traj[0]
    protein.atom_slice(protein.top.select(selection), inplace=True)

    residue_positions = pd.Series([str(x)[3:] for x in list(protein.top.residues)])

    if not os.path.exists("xmgrace"):
        os.makedirs("xmgrace")

    values = pd.Series(val.flatten())
    df = pd.DataFrame({"residues":residue_positions, "values":values})
    df.to_csv("xmgrace/{}_avg_{}.dat".format(prefix, output), header=False, sep="\t", index=False)
