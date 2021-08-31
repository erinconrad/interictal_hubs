function spikes = add_locs_hup132(spikes)

folder = '../../../../../implant_effect/data/elec_locs/HUP132/';
locs_file = [folder,'electrodenames_coordinates_mni.csv'];
ana_file = [folder,'anatomy.csv'];

% Load locs file and ana file
Tl = readtable(locs_file);
Ta = readtable(ana_file,'ReadVariableNames',false);

loc_name = clean_labels_2(Tl.Var1);
locs = [Tl.Var2,Tl.Var3,Tl.Var4];

ana_name = clean_labels_2(Ta.Var1);
ana = Ta.Var2;

for f = 1:length(spikes.file)
    labels = clean_labels_2(spikes.file(f).block(1).chLabels);
    
    % initialize the locs and ana for htat file
    curr_locs = nan(length(labels),3);
    curr_ana = cell(length(labels),1);
    
    % find the indices in loc_name that correspond to labels
    [lia,locb] = ismember(labels,loc_name);
    
    % get the locations corresponding
    curr_locs(lia == 1,:) = locs(locb(lia == 1),:);
    
    % attempt to reconstruct labels from loc_name
    if ~isequal(labels(lia==1),loc_name(locb(lia==1))), error('why'); end
    
 
    % find the indices in ana_name that correpsond to labels
    [lia,locb] = ismember(labels,ana_name);
    curr_ana(lia==1) = ana(locb(lia==1));
    if ~isequal(labels(lia==1),ana_name(locb(lia==1))), error('why'); end
    
    % add to spikes
    spikes.file(f).locs = curr_locs;
    spikes.file(f).ana = curr_ana;
    
end

end