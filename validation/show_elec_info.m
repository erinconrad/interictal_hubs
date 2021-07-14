function show_elec_info(out)

for i = 1:length(out)
    
    unchanged_labels = out(i).unchanged_labels;
    added_labels = out(i).added_labels;
    
    ekg_unchanged = identify_ekg_scalp(unchanged_labels);
    ekg_added = identify_ekg_scalp(added_labels);
    
    unchanged_labels(ekg_unchanged) = [];
    added_labels(ekg_added) = [];
    
    fprintf('\n\n\n%s unchanged:\n',out(i).name);
    unchanged_labels
    
    fprintf('\n\n\n%s added:\n',out(i).name);
    added_labels
    
    
    
end

end