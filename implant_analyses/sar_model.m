function [coeffs,pvals] = sar_model(x,y,locs,pc,which)

wij = get_weights(locs,pc,which);

f = isnan(x) | isnan(y) | isnan(nanmean(wij,2));
x(f) = [];
y(f) = [];
wij(f,:) = [];
wij(:,f) = [];

if isempty(y)
    coeffs = [nan;nan];
    pvals = [nan; nan];
    return
end

try
    results = sar(y,x,wij);
    coeffs = [results.beta;results.rho];
    pvals = norm_prb(results.tstat);
catch
    coeffs = [nan;nan];
    pvals = [nan;nan];
    return
end

%prt(results)

end


function wij = get_weights(locs,pc,which)

if strcmp(which,'pc')
    wij = pc;
    for i = 1:size(wij,1)
        wij(i,i) = 0;
    end
else
    nChannels = size(locs,1);
    dmin = ceil(nanmedian(vecnorm(diff(locs,1),2,2)));

    % Calculate dij, distances between channels and weights wij
    dij = nan(nChannels,nChannels);
    wij = nan(nChannels,nChannels);
    for iChannel = 1:nChannels
       for jChannel = 1:nChannels
           if iChannel == jChannel
               dij(iChannel,jChannel) = 0;
               wij(iChannel,jChannel) = 0;
           else
               % calculate Euclidean distance between the two channels
               dij(iChannel,jChannel) =  vecnorm(locs(iChannel,:)-locs(jChannel,:));

               % if this distance is less than a minimum distance
               if dij(iChannel,jChannel) <= dmin

                   if strcmp(which,'dist1') == 1

                       wij(iChannel,jChannel) = 1;

                   elseif strcmp(which,'dist2') == 1

                       % make the weight be 1/d
                       wij(iChannel,jChannel) = 1/dij(iChannel,jChannel);

                       if wij(iChannel,jChannel) == inf
                           wij(iChannel,jChannel) = 0;
                           fprintf('Warning, channels %d and %d have the same location, making spatial weight 0\n',iChannel,jChannel);
                       end

                   end

               else

                   % if they're further, set the weight to 0
                   wij(iChannel,jChannel) = 0;
               end
           end

       end
    end
    
end

end