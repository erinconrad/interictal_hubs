import subprocess
import sys

def run_matlab(whichPts):

    if len(whichPts) != 2:
        print("Need two patients exactly")
        raise

    txt = '"matlab","-nodisplay","-nodesktop","-r",\
                             "addpath(genpath(''/mnt/local/gdrive/public/USERS/erinconr/projects/interictal_hubs/tools/''));'\
    'cd ../get_hubs;dbstop if error;get_spikes([{} {}]);exit"'.format(whichPts[0],whichPts[1])
    
    subprocess.call([txt])

if __name__ == '__main__':
    whichPts = sys.argv
    run_matlab(whichPts)
