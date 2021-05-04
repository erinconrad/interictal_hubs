function [cos,unchanged_spikey_labels] = co_spike_index(rate_post,spikey_idx,coa_post,blocks,added_labels,unchanged_labels,post_labels,new_post_labels)

%% Find unchanged electrodes that are spikey enough
unchanged_spikey_labels = unchanged_labels(spikey_idx);

%% Get indices of unchanged (that are spikey enough) and added elecs
keep = ismember(post_labels,unchanged_spikey_labels) | ismember(post_labels,added_labels); % only keep spikey unchanged or post

%% Spikey post-labels - remove non-spikey unchanged from post labels
spikey_post_labels = post_labels(keep);
new_post_labels = new_post_labels(keep);
unchanged = find(ismember(spikey_post_labels,unchanged_spikey_labels));
added = find(ismember(spikey_post_labels,added_labels));

%% restrict coa and rate_unchanged to time blocks of interest and electrodes of interest
coa_post = coa_post(keep,keep,blocks); % just spikey-unchanged or post
rate_unchanged = rate_post(spikey_idx,blocks); % just spikey unchanged; note rate_post only shows the rate for unchanged

%% For each unchanged electrode, identify its maximally co-occurring added electrode and how much they co-spike
nunchanged = length(unchanged);
all_max_add = nan(nunchanged,1);
all_max_coa = nan(nunchanged,1);
sum_coa = nan(nunchanged,1);
for i = 1:nunchanged
    ich = unchanged(i);
    
    % co-activation array between ich and all added electrode
    coa_added = coa_post(ich,added,:);
    
    % Take mean across blocks
    coa_added = nanmean(coa_added,3);

    % this is the index of the added electrode with max co-occurrence for
    % this unchanged electrode
    [max_coa,max_add] = max(coa_added);
    max_add = added(max_add);
    all_max_add(i) = max_add;
    all_max_coa(i) = max_coa;
    
    % Sum across all other electrodes to get the sum coa
    sum_coa(i) = nansum(nanmean(coa_post(ich,:,:),3));
end

%% Divide this by its overall co-spikeyness to get its relative co-spike index
%cos = all_max_coa./nanmean(rate_unchanged,2);
cos = all_max_coa./sum_coa;

cos_all = nanmean(coa_post(:,:,blocks),3);
cos_all_rel = cos_all./repmat(sum_coa,1,size(cos_all,2));

%% show stuff
if 1
    figure
    set(gcf,'position',[0 0 1300 800]);
    %imagesc(nanmean(coa_post(:,:,blocks),3))
    imagesc(cos_all)
    xticks(1:length(spikey_post_labels))
    yticks(1:length(spikey_post_labels))
    xticklabels(new_post_labels)
    yticklabels(new_post_labels)
    xtickangle(90)
    
end

end