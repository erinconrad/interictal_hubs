function [median_rl,nseq] = get_sequences(gdf,nchs)

t2 = 15*1e-3; % max time from preceding spike (15 ms in paper)
minSeqLength = 5; 
nspikes = size(gdf,1);

if nspikes == 0
    mean_rl = nan(nchs,1);
    median_rl = nan(nchs,1);
    nseq = 0;
    return
end
 
%% Confirm they're ordered in time
times = gdf(:,2);
if ~isequal(times,sort(times))
    error('ruh-roh');
end

%% initialize seq cell
seq = {};

%% Initialize stuff for my while loop
curr_seq = gdf(1,:);
next_index = 2;

while next_index <= nspikes
    
    next_time = gdf(next_index,2);
    last_time = curr_seq(end,2);
    
    % go to the next if the channel is already in the sequence
    if ismember(gdf(next_index,1),curr_seq(:,1))
        next_index = next_index+1;
        continue;
    end
    
    % if they're close enough in time
    if next_time - last_time < t2
        % append it to the sequence
        curr_seq = [curr_seq;gdf(next_index,:)];
        
        % increase next_index
        next_index = next_index + 1;
        
    % if the next one is too far
    else
        
        
        % check if the current sequence is long enough
        if size(curr_seq,1) >= minSeqLength
            % Append it to the seq cell array
            seq = [seq;curr_seq];
            
            % Reset curr seq to be this new spike
            curr_seq = gdf(next_index,:);
            
            % Increase next index
            next_index = next_index+1;
            
        % If it's too short, don't add it but increase index
        else
            % Reset curr seq to be this new spike
            curr_seq = gdf(next_index,:);
            
            % Increase next index
            next_index = next_index+1;
            
        end
        
    end
    
end

nseq = length(seq);

%% Now that I have my sequences, get some info
rl = nan(nchs,nseq);
degree = nan(nchs,nseq); % note this is probably crap because of how many are detected simultaneously

for s = 1:nseq
    for ich = 1:nchs
        curr = seq{s};
        
        if ismember(ich,curr(:,1))
            
            ich_index = find(curr(:,1) == ich);
            
            % difference between the time this electrode peaks and the
            % first peak in the sequence
            rl(ich,s) = curr(ich_index,2) - curr(1,2);
            
            % number after - number before
            degree(ich,s) =  (size(curr,1) - ich_index) - (ich_index - 1);
        end
    end
end

mean_rl = nanmean(rl,2);
median_rl = nanmedian(rl,2);
meancorr = nan(nseq,1);
mediancorr = nan(nseq,1);

if 0
    for s = 1:nseq
        meancorr(s) = corr(mean_rl,rl(:,s),'Type','Spearman','rows','pairwise');
        mediancorr(s) = corr(median_rl,rl(:,s),'Type','Spearman','rows','pairwise');
    end
    
    plot(meancorr,'o')
    hold on
    plot(mediancorr,'x')
end

end