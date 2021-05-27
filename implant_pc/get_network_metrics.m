function out = get_network_metrics(net)

%% Derive sizes
nblocks = size(net,2);
nchs = (1+sqrt(1+8*size(net,1)))/2;


%% initialize
ns_norm = nan(nchs,nblocks);
ge = nan(1,nblocks);

%% Loop over blocks
for ib = 1:nblocks
    
    %% unpack matrix
    mat = wrap_or_unwrap_adjacency(net(:,ib));
    
    %% Get network metrics
    %{
    Note that because some channels are nans and which are nans change at
    different times, I need to normalize!
    %}
    
    % Node strength normalized
    ns_norm_temp = nanmean(mat,1);
    ns_norm(:,ib) = ns_norm_temp;
    
    % Global efficiency - I believe this is already normalized, need to
    % check!
    E = efficiency_wei(mat,0);
    ge(ib) = E;
    
end

out.ns_norm = ns_norm;
out.ge = ge;

end