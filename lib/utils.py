import os

import numpy as np
import mdtraj as md



class MDIterator:

    def __init__(self, traj_file, top, chunk=100, stride=1):
        self.iterator = md.iterload(traj_file, top=top, chunk=chunk, stride=stride)
        self.trajectory = None

        self.index = chunk - 1
        self.chunk = chunk


    def __iter__(self):
        return self


    def next(self):
        if self.index == self.chunk - 1:
            self.index = -1
            self.trajectory = self.iterator.next()

        self.index += 1

        try:
            return self.trajectory[self.index]
        except IndexError as ex:
            raise StopIteration



def reduce_trajectory(trajectory, top=None, stride=1, output_path="minimized.dcd"):
    traj = md.load(trajectory, top=top)[::int(stride)]
    traj.save(output_path)


def save_frame(trajectory, topology, frame_index, format="pdb"):
    traj_wo_ext = ".".join(os.path.basename(trajectory).split(".")[:-1])
    frame_name = "%s_%d.%s" % (traj_wo_ext, frame_index, format)
    frame = md.load_frame(trajectory, frame_index, top=topology)
    frame.save(frame_name)


def format_seconds(seconds):
    m, s = divmod(seconds, 60)
    h, m = divmod(m, 60)
    return "%d:%02d:%02d" % (h, m, s)


def saveNDtxt(filename, data):
    with file(filename, 'w') as outfile:
        # I'm writing a header here just for the sake of readability
        # Any line starting with "#" will be ignored by numpy.loadtxt
        outfile.write('# Array shape: {0}\n'.format(data.shape))

        # Iterating through a ndimensional array produces slices along
        # the last axis. This is equivalent to data[i,:,:] in this case
        for data_slice in data:

            # The formatting string indicates that I'm writing out
            # the values in left-justified columns 7 characters in width
            # with 2 decimal places.
            np.savetxt(outfile, data_slice, fmt='%-7.2f')

            # Writing out a break to indicate different slices...
            outfile.write('# New slice\n')


def loadNDtxt(filename, shape):
    return np.loadtxt(filename).reshape(shape)

def dat2xmgrace(val, prefix, output, traj, selection="(name CB and protein) or (name CA and resname GLY)"):
    import pandas as pd
    assert type(traj) is md.core.trajectory.Trajectory, "traj has to be an object of type md.core.trajectory.Trajectory"

    protein = traj[0]
    protein.atom_slice(protein.top.select(selection), inplace=True)

    residue_positions = pd.Series([str(x)[3:] for x in list(protein.top.residues)])

    import os
    if not os.path.exists("xmgrace"):
        os.makedirs("xmgrace")

    values = pd.Series(val.flatten())
    df = pd.DataFrame({"residues":residue_positions, "values":values})
    df.to_csv("xmgrace/{}_avg_{}.dat".format(prefix, output), header=False, sep="\t", index=False)

'''
Example usage:
'''
if __name__ == "__main__":
    import sys
    '''
    traj = MDIterator(sys.argv[1], sys.argv[2])
    print dir(traj)

    for i, frame in enumerate(traj):
        print i
    '''

    reduce_trajectory(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])




