clear

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];
ieeg_folder = locations.ieeg_folder;
pwfile = locations.ieeg_pw_file;
login_name = locations.ieeg_login;

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

% Loop over patients
for p =1:length(pt)
    
    name = pt(p).name;
    fprintf('\nDoing %s\n',name);
    
    if strcmp(name,'HUP060')
        base_ieeg_name = [name,'_phaseIV'];
    else
        base_ieeg_name = [name,'_phaseII'];
    end
    
    dcount = 0;
    add_it = 0;
    finished = 0;
    
    while 1
        if dcount == 0
            ieeg_name = base_ieeg_name;
            try
                session = IEEGSession(ieeg_name,login_name,pwfile);
                finished = 1;
                add_it = 1;
                dcount = 1;
            catch
                try
                    % First, remove leading 0
                    [s,e] = regexp(name,'\d*');
                    num_str = name(s:e);
                    if strcmp(num_str(1),'0')
                        new_name = [name(1:s-1),name(s+1:e)];
                        ieeg_name = [new_name,'_phaseII'];
                        session = IEEGSession(ieeg_name,login_name,pwfile);
                        finished = 1;
                        add_it = 1;
                        dcount = 1;
                    else
                        error('')
                    end
                    
                catch
                    fprintf('\nDid not find %s, adding an appendage\n',ieeg_name);
                end
            end
        else
            
            ieeg_name = [base_ieeg_name,'_D0',sprintf('%d',dcount)];
            try
                session = IEEGSession(ieeg_name,login_name,pwfile);
                finished = 0;
                add_it = 1;
            catch
                add_it = 0;
                finished = 1; % if I can't find it adding
            end
        end
        
        
        % Add session info
        if add_it == 1
            pt(p).ieeg.file(dcount).fs = session.data.sampleRate;
            pt(p).ieeg.file(dcount).name = session.data.snapName;
            pt(p).ieeg.file(dcount).chLabels = session.data.channelLabels;
            pt(p).ieeg.file(dcount).duration = session.data.rawChannels(1).get_tsdetails.getDuration;
        end
        
        if finished == 1
            if exist('session','var') ~= 0
                session.delete;
            end
            break
        end
        
        dcount = dcount + 1;
        if exist('session','var') ~= 0
            session.delete;
        end
    end
    
end

% Save the pt struct
save([data_folder,'pt.mat'],'pt');
