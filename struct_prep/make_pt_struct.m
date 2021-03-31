clear

%% Get file locs
locations = interictal_hub_locations;
data_folder = [locations.main_folder,'data/'];

%% Load json file into struct
data = loadjson([data_folder,'DATA_MASTER.json']);

ptnames = fieldnames(data.PATIENTS);

for p = 1:length(ptnames)
    
    curr = data.PATIENTS.(ptnames{p});
    
    % add name
    pt(p).name = ptnames{p};
    
    % add clinical info
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
        pt(p).seizure_info.sz(i).EEC = curr_sz.SeizureEEC;
        pt(p).seizure_info.sz(i).UEO = curr_sz.SeizureUEO;
        if isfield(curr_sz,'SeizureEnd')
            pt(p).seizure_info.sz(i).End = curr_sz.SeizureEnd;
        else
            pt(p).seizure_info.sz(i).End = [];
        end
        pt(p).seizure_info.sz(i).OnsetElectrodes = curr_sz.SEIZURE_ONSET_ELECTRODES;
        
    end
    
    
end