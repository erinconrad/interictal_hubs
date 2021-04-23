function res_idx = resected_ch_nums(resected,labels)

%% Clean both sets
clean_labels = clean_labels_2(labels);
clean_resected = clean_labels_2(resected);

%% Get channel indices of resected electrodes
res_idx = find(ismember(clean_labels,clean_resected));

end