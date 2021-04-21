function [seq,rl] = new_get_sequences(gdf,nchs)
    
t2 = 15*1e-3; % max time from preceding spike (15 ms in paper)
minSeqLength = 5; 


ns = size(gdf,1);
spike = 1:ns;

%% Sort by time
times = gdf(:,2);
if ~isequal(times,sort(times))
    error('ruh-roh');
end

%% initialize seq cell
seq = {};

%% Initialize stuff for my while loop
curr_seq = [gdf(1,:),1];
next_index = 2;

while next_index <= ns
    
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
        curr_seq = [curr_seq;[gdf(next_index,:) next_index]];
        
        % increase next_index
        next_index = next_index + 1;
        
    % if the next one is too far
    else
        
        
        % check if the current sequence is long enough
        if size(curr_seq,1) >= minSeqLength
            % Append it to the seq cell array
            seq = [seq;curr_seq];
            
            % Reset curr seq to be this new spike
            curr_seq = [gdf(next_index,:) next_index];
            
            % Increase next index
            next_index = next_index+1;
            
        % If it's too short, don't add it but increase index
        else
            % Reset curr seq to be this new spike
            curr_seq = [gdf(next_index,:) next_index];
            
            % Increase next index
            next_index = next_index+1;
            
        end
        
    end
    
end

nseq = length(seq);

%% Now that I have my sequences, get some info
rl = nan(nchs,nseq);

for s = 1:nseq
    curr = seq{s};
    for ich = 1:nchs
        
        
        if ismember(ich,curr(:,1))
            
            ich_index = find(curr(:,1) == ich);
            
            % difference between the time this electrode peaks and the
            % first peak in the sequence
            rl(ich,s) = curr(ich_index,2) - curr(1,2);
            
        end
    end
end



end