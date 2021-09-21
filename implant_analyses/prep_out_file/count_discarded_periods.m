function count_discarded_periods(out)

all_perc_discarded = nan(length(out),1);

for p = 1:length(out)
    
    rate = out(p).rate;
    
    % get number of periods where all electrodes have nan
    discarded_periods = sum(isnan(rate),1) == size(rate,1);
    
    perc_discarded = sum(discarded_periods)/length(discarded_periods);
    all_perc_discarded(p) = perc_discarded;
    
end

fprintf(['\nAn average of %1.1f%% segments (range %1.1f%%-%1.1f%%) total were rejected as either '...
    'apparently disconnected or artifact-heavy by this method.\n'],...
    mean(all_perc_discarded)*100,min(all_perc_discarded)*100,max(all_perc_discarded)*100);

end