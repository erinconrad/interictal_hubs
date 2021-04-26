function new_spike_collection(whichPts)

if length(whichPts) ~= 2, error('Enter two patients'); end

str = sprintf(['screen & matlab -nodisplay -r '...
    '"addpath(genpath(''/mnt/local/gdrive/public/USERS/erinconr/projects/interictal_hubs/tools/''));'...
    'cd ../get_hubs;get_spikes([%d %d]);exit"'],whichPts(1),whichPts(2));

system(str);

end