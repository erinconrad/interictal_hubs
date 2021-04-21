

%{
This is a one time use function to add patients from the seizure spread
project to the pt structure

%}

clear

%% Get file locs
locations = interictal_hub_locations;
pt_folder = [locations.script_folder,'data/'];
json_folder = [locations.main_folder,'data/'];


%% Load current pt structure
pt = load([pt_folder,'pt.mat']);
pt = pt.pt;

npts = length(pt);
p = npts+1;

%% Load json
data = loadjson([json_folder,'iEEGdataRevell.json']);

% Loop over patients in json
ptnames = fieldnames(data.SUBJECTS);

for q = 1:length(ptnames)
    curr = data.SUBJECTS.(ptnames{q});
    hupname = curr.HUP;
    
    %% Skip if the hup name is missing
    if strcmp(hupname,'missing')
        continue;
    end
    
    %% Skip if not a number
    if ~strcmp(class(hupname),'double')
        continue;
    end
    
    %% turn the hup number into a name
    if hupname < 100
        hupstr = sprintf('HUP0%d',hupname);
    else
        hupstr = sprintf('HUP%d',hupname);
    end
    
    %% Skip if I already have this patient
    found_it = 0;
    for j = 1:length(pt)
        if strcmp(hupstr,pt(j).name)
            found_it = 1;
            break
        end
    end
    if found_it == 1
        continue;
    end
    
    %% If I don't have it yet, add it to struct
    
    % Add name
    pt(end+1).name = hupstr;
    
    % Add clinical info
    pt(p).clinical.sex = curr.Sex;
    
    if isfield(curr,'Outcome_0x20_6_0x20_month')
        pt(p).clinical.outcomes.month6 = curr.Outcome_0x20_6_0x20_month;
    else
        pt(p).clinical.outcomes.month6 = '';
    end
    
    if isfield(curr,'Outcome_0x20_12_0x20_month')
        pt(p).clinical.outcomes.month12 = curr.Outcome_0x20_12_0x20_month;
    else
        pt(p).clinical.outcomes.month12 = '';
    end
    
    if isfield(curr,'Outcome_0x20_24_0x20_month')
        pt(p).clinical.outcomes.month24 = curr.Outcome_0x20_24_0x20_month;
    else
        pt(p).clinical.outcomes.month24 = '';
    end
    
    pt(p).clinical.AgeOnset = curr.AgeOnset;
    pt(p).clinical.AgeSurgery = curr.AgeSurgery;
    
    if isfield(curr,'Location')
        pt(p).clinical.Location = curr.Location;
    else
        pt(p).clinical.Location = '';
    end
    pt(p).clinical.LesionStatus = curr.Lesion_0x20_Status;
    
    if isfield(curr,'Previous_0x20_Surgery')
        pt(p).clinical.PreviousSurgery = curr.Previous_0x20_Surgery;
    else
        pt(p).clinical.PreviousSurgery = '';
    end
    pt(p).clinical.Pathology = curr.Pathology;
    
    if isfield(curr,'Surgery_0x20_Type')
        pt(p).clinical.SurgeryType = curr.Surgery_0x20_Type;
    else
        pt(p).clinical.SurgeryType = '';
    end
    
    % add some electrode info
    if isfield(curr,'RESECTED_ELECTRODES')
        pt(p).electrode_info.resected = curr.RESECTED_ELECTRODES;
    elseif isfield(curr,'MANUAL_RESECTED_ELECTRODES')
        pt(p).electrode_info.resected = curr.MANUAL_RESECTED_ELECTRODES;
    else
        pt(p).electrode_info.resected = {};
    end
    pt(p).electrode_info.ignore = curr.IGNORE_ELECTRODES;
    
    % Add events info
    sznames = fieldnames(curr.Events.Ictal);
    for i = 1:length(sznames)
        curr_sz = curr.Events.Ictal.(sznames{i});
        pt(p).seizure_info.sz(i).type = curr_sz.SeizureType;
        pt(p).seizure_info.sz(i).EMU_event_no = curr_sz.EMU_Report_Event_Number;
        if isfield(curr_sz,'EEC')
            pt(p).seizure_info.sz(i).EEC = curr_sz.EEC;
        else
            pt(p).seizure_info.sz(i).EEC = [];
        end
        pt(p).seizure_info.sz(i).UEO = curr_sz.UEO;
        if isfield(curr_sz,'Stop')
            pt(p).seizure_info.sz(i).End = curr_sz.Stop;
        else
            pt(p).seizure_info.sz(i).End = [];
        end
        pt(p).seizure_info.sz(i).OnsetElectrodes = curr_sz.OnsetElectrodes;
        
    end
    
    %% Add one to the latest index
    p = p+1;
    
end

