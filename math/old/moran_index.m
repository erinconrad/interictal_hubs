function MI = moran_index(locs,data,dmin)
data = data';
% set dmin

wij = getwij2(locs,dmin);

N = length(data);
E = -1/(N-1);


xdmean = data-nanmean(data);
W = nansum(wij(:));



op = xdmean' * xdmean;
if sum(size(op)) == 2
    error('Warning, I think you need to take the transpose of rl/n');
end

I = N/W*nansum(nansum(wij.*op))/...
    nansum((xdmean).^2);
S1 = 1/2*nansum(nansum((wij+wij').^2));
S2 = nansum((nansum(wij)+nansum(wij')).^2);
S3 = 1/N*nansum((xdmean).^4)/(1/N*nansum((xdmean).^2))^2;
S4 = (N^2-3*N+3)*S1 - N*S2 + 3*W^2;
S5 = (N^2-N)*S1 - 2*N*S2 + 6*W^2;

V = (N*S4-S3*S5)/((N-1)*(N-2)*(N-3)*W^2)-E^2;

Z = (I-E)/sqrt(V);
p = normcdf(-Z);

MI.I = I;
MI.E = E;
MI.V = V;
MI.Z = Z;
MI.p = p;
MI.wij = wij;

end

function wij = getwij2(locs,dmin)

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
           dij(iChannel,jChannel) =  sqrt(sum((locs(iChannel,:) - locs(jChannel,:)).^2));
           
           % if this distance is less than a minimum distance
           if dij(iChannel,jChannel) <= dmin
               
               % make the weight be 1/d
               wij(iChannel,jChannel) = 1/dij(iChannel,jChannel);
               
               if wij(iChannel,jChannel) == inf
                   wij(iChannel,jChannel) = 0;
                   fprintf('Warning, channels %d and %d have the same location, making spatial weight 0\n',iChannel,jChannel);
               end
           else
               
               % if they're further, set the weight to 0
               wij(iChannel,jChannel) = 0;
           end
       end
      
   end
    
end




end
