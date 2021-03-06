function out = get_network_metrics(net,unchanged_labels,rm_ekg)

ekg = identify_ekg_scalp(unchanged_labels);

%% Derive sizes
nblocks = size(net,2);
nchs = (1+sqrt(1+8*size(net,1)))/2;


%% initialize
ns_norm = nan(nchs,nblocks);
ns = nan(nchs,nblocks);
ge = nan(1,nblocks);
all_mat = nan(nchs,nchs,nblocks);

%% Loop over blocks
for ib = 1:nblocks
    
    %% unpack matrix
    mat = wrap_or_unwrap_adjacency(net(:,ib));
    mat(1:1+size(mat,1):end) = nan;
    
    %% Get network metrics
    %{
    Note that because some channels are nans and which are nans change at
    different times, I need to normalize!
    %}
    % set EKG to nans
    if rm_ekg
        mat(ekg,:) = nan;
        mat(:,ekg) = nan;
    end
        
    
    % Node strength normalized
    ns_norm_temp = nanmean(mat,1);
    ns_norm(:,ib) = ns_norm_temp;
    
    % Absolute node strength
    ns_temp = nansum(mat,1);
    ns_temp(isnan(nanmean(mat,1))) = nan; % turn things that are all nan into nan
    ns(:,ib) = ns_temp;
    
    % Global efficiency - I believe this is already normalized, need to
    % check!
    E = efficiency_wei(mat,0);
    ge(ib) = E;
    
    % all mat
    all_mat(:,:,ib) = mat;
    
end

out.ns_norm = ns_norm;
out.ns = ns;
out.ge = ge;
out.avg_mat = nanmean(all_mat,3);

end