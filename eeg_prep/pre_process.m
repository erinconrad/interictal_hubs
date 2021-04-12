function values = pre_process(values,chLabels)

nchs = size(values,2);

%% Bipolar montage
for ch = 1:nchs
    
    % Initialize it as nans
    out = nan(size(values,1),1);

    % Get electrode name
    label = chLabels{ch};

    % get the non numerical portion
    label_num_idx = regexp(label,'\d');
    label_non_num = label(1:label_num_idx-1);

    % get numerical portion
    label_num = str2num(label(label_num_idx:end));

    % see if there exists one higher
    label_num_higher = label_num + 1;
    higher_label = [label_non_num,sprintf('%d',label_num_higher)];
    if sum(strcmp(chLabels(:,1),higher_label)) > 0
        higher_ch = find(strcmp(chLabels(:,1),higher_label));
        out = values(:,ch)-values(:,higher_ch);
    else
        % allow it to remain nans
    end
    values(:,ch) = out;
    
    
end


end