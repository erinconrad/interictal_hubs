function [which_chs,skip] = designate_chs(chLabels)

which_chs = 1:length(chLabels);
exclusion_text = {'ekg','EKG','ecg','ECG','rate','rr','RR'};

to_exclude = cellfun(@(x) contains(x,exclusion_text),chLabels,'UniformOutput',false);
to_exclude = cell2mat(to_exclude);

skip = which_chs(to_exclude);

which_chs(to_exclude) = [];


end