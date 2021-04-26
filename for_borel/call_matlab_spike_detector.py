import subprocess
import sys

def run_matlab(whichPts):
    print(whichPts)
    print(type(whichPts))
    if len(whichPts) != 2:
        print("Need two patients exactly")
        raise

    txt = '"matlab","-nodisplay","-nodesktop","-r",\
                             "dbstop if error;fprintf(''Hello'');addpath(genpath(''/mnt/local/gdrive/public/USERS/erinconr/projects/interictal_hubs/tools/''));'\
    'cd ../get_hubs;get_spikes([{} {}]);exit"'.format(whichPts[0],whichPts[1])
    
    subprocess.call([txt])

if __name__ == '__main__':
    whichPts = [sys.argv[1],sys.argv[2]]
    run_matlab(whichPts)
