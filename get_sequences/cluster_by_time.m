function thing_clusts = cluster_by_time(thing,time,cluster_time)

clusts = time(1):cluster_time:time(end);
nclusts = length(clusts);
thing_clusts = nan(size(thing,1),nclusts);

for i = 1:nclusts-1
    curr_times = [clusts(i) clusts(i+1)];
    % find indices in this time
    in_clust = time >= curr_times(1) & time < curr_times(2);
    
    clust_avg = nanmean(thing(:,in_clust),2);
    thing_clusts(:,i) = clust_avg;

end

end