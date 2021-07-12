function elec_details(out)

num_added = nan(length(out),1);
num_unchanged = nan(length(out),1);
names = cell(length(out),1);

for i = 1:length(out)
    added_labels = out(i).added_labels;
    ekg_added = identify_ekg_scalp(added_labels);
    added_labels(ekg_added) = [];
    
    unchanged_labels = out(i).unchanged_labels;
    ekg_unchanged = identify_ekg_scalp(unchanged_labels);
    unchanged_labels(ekg_unchanged) = [];
    
    name = out(i).name;
    
    if 0
        fprintf('\n%s unchanged:\n',out(i).name);
        unchanged_labels
        
        fprintf('\n%s added:\n',out(i).name);
        added_labels
        
        pause
    end
    
    num_added(i) = length(added_labels);
    num_unchanged(i) = length(unchanged_labels);
    names{i} = name;
    
end

fprintf(['\nThere was an average of %1.1f original electrodes (range %d-%d), '...
    'and %1.1f added electrodes (range %d-%d).\n'],...
    mean(num_unchanged),min(num_unchanged),max(num_unchanged),...
    mean(num_added),min(num_added),max(num_added));

num_unchanged

num_added
end