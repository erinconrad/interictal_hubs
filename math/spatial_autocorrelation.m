function [coeffs,pvals] = spatial_autocorrelation(x,y,dmin,locs,un_cos,first_order)
wij = getwij3(locs,dmin,un_cos,first_order);

%{
A = zeros(length(y),1);
for i = 1:size(A,1)
    for j = 1:size(A,1)
        A(i) = A(i) + wij(i,j)*y(j);
    end
end

tbl = table(y,x,A);
lm = fitlm(tbl,'y~x+A');
%}

% Remove nans
f = isnan(x) | isnan(y);
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
    %}

%prt(results)

end


function wij = getwij3(locs,dmin,un_cos,first_order)


if first_order == 2
    
    wij = un_cos;
    
else

    nChannels = size(locs,1);

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

                   if first_order == 1

                       wij(iChannel,jChannel) = 1;

                   else

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