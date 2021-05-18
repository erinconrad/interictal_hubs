function histogram_rate_change(increase,labels)


%% Remove ekg
ekg = identify_ekg_scalp(labels);
labels(ekg) = [];
increase(ekg) = [];

%% Show histogram of increase
figure
histogram(increase)




end