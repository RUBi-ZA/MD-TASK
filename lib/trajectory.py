import os, math
import mdtraj as md

class MDIterator(object):

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
        except IndexError:
            raise StopIteration

def reduce_trajectory(trajectory, top=None, stride=1, output_path="minimized.dcd"):
    traj = md.load(trajectory, top=top)[::int(stride)]
    traj.save(output_path)

def save_frame(trajectory, topology, frame_index, file_format="pdb"):
    traj_wo_ext = ".".join(os.path.basename(trajectory).split(".")[:-1])
    frame_name = "%s_%d.%s" % (traj_wo_ext, frame_index, file_format)
    frame = md.load_frame(trajectory, frame_index, top=topology)
    frame.save(frame_name)

def load_trajectory(trajectory, topology, step=1, lazy_load=False):
    if not lazy_load:
        traj = md.load(trajectory, top=topology)[::step]
        total_frames = len(traj)
    else:
        traj = MDIterator(trajectory, top=topology, stride=step)
        total_frames = None

    return traj, total_frames

def calc_distance(frame, index1, index2):
    atom1 = frame.xyz[0, index1]
    atom2 = frame.xyz[0, index2]

    dist = math.sqrt((atom2[0] - atom1[0])**2 + (atom2[1] - atom1[1])**2 + (atom2[2] - atom1[2])**2)

    return abs(dist)
