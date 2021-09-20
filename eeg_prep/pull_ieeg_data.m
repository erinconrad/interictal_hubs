function values = pull_ieeg_data(fname, login_name, pwfile, run_idx)

attempt = 1;

% Wrap data pulling attempts in a while loop
while 1
    
    try

        session = IEEGSession(fname, login_name, pwfile);
        channelLabels = session.data.channelLabels;
        nchs = size(channelLabels,1);

        % Break the number of channels in half to avoid wacky server errors
        values1 = session.data.getvalues(run_idx,1:floor(nchs/4));
        values2 = session.data.getvalues(run_idx,floor(nchs/4)+1:floor(2*nchs/4));
        values3 = session.data.getvalues(run_idx,floor(2*nchs/4)+1:floor(3*nchs/4)); 
        values4 = session.data.getvalues(run_idx,floor(3*nchs/4)+1:nchs); 

        

        values = [values1,values2,values3,values4];
        
        % break out of while loop if I got data
        break
    
    % If server error, try again (this is because there are frequent random
    % server errors).
    catch ME
        if contains(ME.message,'203') || contains(ME.message,'204') || ...
                contains(ME.message,'202')
            attempt = attempt + 1;
            fprintf('Failed to retrieve ieeg.org data, trying again (attempt %d)\n',attempt); 
        else
            error('Non-server error');
        end
    end
        
session.delete;
clearvars -except values

end