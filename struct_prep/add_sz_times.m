clear

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.script_folder,'data/'];

%% Load pt struct
pt = load([data_folder,'pt.mat']);
pt = pt.pt;

%% Get sheet names from sz table
sheet = sheetnames([data_folder,'seizure times.xlsx']);

% Loop over patients
for p = 1:length(pt)
    name = pt(p).name;
    
    for s = 1:length(sheet)
        
        sname = sheet(s);
        if contains(sname,name)
            if ~contains(sname,'complete')
                break
            else
                
                % Load sheet
                T = readtable([data_folder,'seizure times.xlsx'],'Sheet',sname);
                
                % initialize data
                for f = 1:length(pt(p).ieeg.file)
                    pt(p).ieeg.file(f).sz = [];
                end
                
                % fill sz info
                for sz = 1:size(T,1)
                    f = T.IEEG_file(sz);
                    info.eec = T.EEC(sz);
                    info.ueo = T.UEO(sz);
                    info.end = T.End(sz);
                    info.elec = T.OnsetElectrode(sz);
                    info.type = T.Type(sz);
                    info.notes = T.Notes(sz);
                    
                    pt(p).ieeg.file(f).sz = [pt(p).ieeg.file(f).sz;
                        info];
                end
                
                % Quit looking for matching pt
                break
            
            end
            
            
        end


    end
    
end