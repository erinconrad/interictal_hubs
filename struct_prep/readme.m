%{
This folder contains scripts to prep the data structure. It should only
need to be done once, or when making major changes to how I will run the
project. The steps of this include (in order):

1) make_pt_struct: this takes the DATA_MASTER.json file and converts it to
a patient struct
2) identify_ieeg_files: this takes the pt struct and fills up info about
the different ieeg files, fs, durations, and chlabels
3) define_times: this takes the pt struct populated with ieeg info, and
picks random times over which the spike detector will be run. 
    - this calls intersect_sz_time to avoid seizures in defining run times


Miscellaneous scripts:
- clean_up_extra_entries: this was made to fix an early bug and should not
need to be used again

%}