function qual_analysis

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
qual_results = [results_folder,'qual/'];
if ~exist(qual_results,'dir')
    mkdir(qual_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.main_folder,'data/'];
dname = [data_folder,'reimplantation patients.xlsx'];

T = readtable(dname,'sheet','Qual table 2');

%% Get response and predictor variable
surg = contains(T.InterventionOffered,'ATL') | contains(T.InterventionOffered,'Ablation');
capture = contains(T.RevisionCapturedOnset,'Yes');
orig_elec_gs = contains(T.OriginalElectrodes,'S') | contains(T.OriginalElectrodes,'G');
revised_elec_gs = contains(T.RevisedElectrodes,'S') | contains(T.RevisedElectrodes,'G');
semlat = contains(T.LateralizingSemiology,'Yes'); % NEED TO GET FINAL PATIENTS
mriles = contains(T.LesionalMRI,'Yes');
petles = contains(T.LesionalPET,'Yes');
rev_reas = T.ReasonForRevision; % categorical
missed_target = contains(T.ReasonForRevision,'missed');
number = T.Subject;

%% Get stats for binary-binary tests
% Orig electrode type -> surg offered
tbl = crosstab(orig_elec_gs,surg);
[~,p_orig,stats] = fishertest(tbl);
or_orig = stats.OddsRatio;
orig_text = sprintf('Original electrodes\tStereo-EEG (%d/%d)\tOther (%d/%d)\t%1.1f\t%1.2f',...
    sum(orig_elec_gs==0&surg==1),sum(orig_elec_gs==0),...
    sum(orig_elec_gs==1&surg==1),sum(orig_elec_gs==1),...
    or_orig,p_orig);

% Revised electrode type -> surg offered
tbl = crosstab(revised_elec_gs,surg);
[~,p_rev,stats] = fishertest(tbl);
or_rev = stats.OddsRatio;
rev_text = sprintf('Revised electrodes\tStereo-EEG (%d/%d)\tOther (%d/%d)\t%1.1f\t%1.2f',...
    sum(revised_elec_gs==0&surg==1),sum(revised_elec_gs==0),...
    sum(revised_elec_gs==1&surg==1),sum(revised_elec_gs==1),...
    or_rev,p_rev);

% Lateralized semiology -> surg offered
tbl = crosstab(semlat,surg);
[~,p_sem,stats] = fishertest(tbl);
or_sem = stats.OddsRatio;
sem_text = sprintf('Lateralizing semiology\tYes (%d/%d)\tNo (%d/%d)\t%1.1f\t%1.2f',...
    sum(semlat==1&surg==1),sum(semlat==1),...
    sum(semlat==0&surg==1),sum(semlat==0),...
    or_sem,p_sem);

% Lesional MRI -> surg offered
tbl = crosstab(mriles,surg);
[~,p_mri,stats] = fishertest(tbl);
or_mri = stats.OddsRatio;
mri_text = sprintf('Lesional MRI\tYes (%d/%d)\tNo (%d/%d)\t%1.1f\t%1.2f',...
    sum(mriles==1&surg==1),sum(mriles==1),...
    sum(mriles==0&surg==1),sum(mriles==0),...
    or_mri,p_mri);

% Lesional PET -> surg offered
tbl = crosstab(petles,surg);
[~,p_pet,stats] = fishertest(tbl);
or_pet = stats.OddsRatio;
pet_text = sprintf('Lesional PET\tYes (%d/%d)\tNo (%d/%d)\t%1.1f\t%1.2f',...
    sum(petles==1&surg==1),sum(petles==1),...
    sum(petles==0&surg==1),sum(petles==0),...
    or_pet,p_pet);

% Missed target -> surg offered
tbl = crosstab(missed_target,surg);
[~,p_missed,stats] = fishertest(tbl);
or_missed = stats.OddsRatio;
missed_text = sprintf('Intended target missed\tYes (%d/%d)\tNo (%d/%d)\t%1.1f\t%1.2f',...
    sum(missed_target==1&surg==1),sum(missed_target==1),...
    sum(missed_target==0&surg==1),sum(missed_target==0),...
    or_missed,p_missed);

%% Make a table showing this
C = sprintf('%s\n%s\n%s\n%s\n%s\n%s',...
    orig_text,rev_text,sem_text,mri_text,pet_text,missed_text);

%% Stats for categorical->binary test (3x2 table)
% revision reason -> surg offered
rev_tbl = crosstab(rev_reas,surg); 
% Need to put this into R because Matlab does not have functionality for a
% 2x3 fisher exact test


%% Stats for rank -> binary test
% Patient number -> surg offered
p_num = ranksum(number(surg == 1),number(surg == 0));

%

end