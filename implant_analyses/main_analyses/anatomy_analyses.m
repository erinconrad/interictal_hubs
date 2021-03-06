function [output_for_plot,out_loc_names] = anatomy_analyses(whichPts,saved_out,out)

%{
This analysis tests whether the relative change in electrographic features
is larger for different anatomical regions.
%}

%% Parameters
all_surrounds = 12*[0.5,1,2,3,4,5,6,7,8,9,10];
do_norm = 0; % doesn't seem to make much difference so I will keep 0 for simplicity
main_surround = 3;
which_resp = 'rel';%'rel';% 
which_test = 'skillmack';%'friedman';

%% Other info
n_surrounds = length(all_surrounds);
all_metrics = {'spike rate','node strength'};
n_metrics = length(all_metrics);

% divide by 2 to get time in hours but then multiply by 2 to get the full
% pre+post surround
true_times = all_surrounds; 

%% Locations
locations = interictal_hub_locations;
results_folder = [locations.main_folder,'results/'];
main_spike_results = [results_folder,'main_spikes/'];
if ~exist(main_spike_results,'dir')
    mkdir(main_spike_results);
end

addpath(genpath(locations.script_folder));
data_folder = [locations.script_folder,'data/'];
spike_folder = [results_folder,'new_spikes/'];

if isempty(whichPts)
    whichPts = [20 103 106 107 35 109 110 111 94 97];
end

if saved_out == 1
    
    %out = load([main_spike_results,'out.mat']);
    %out = out.out;
    
else
    out = initialize_out_struct(length(whichPts));
    
    %% Get spike details
    fprintf('Getting spike details for pt...\n');
    for i = 1:length(whichPts)
        p = whichPts(i);
        fprintf('%d of %d\n',i,length(whichPts));
        out(i) = get_gdf_details(p);
    end
    save([main_spike_results,'out'],'out');
end

%% Get anatomy loc names
[~,~,loc_names] = anatomy_grouper([]);

%% Prep final table
final_array = cell(length(all_metrics),n_surrounds);
loc_nums = zeros(n_metrics,n_surrounds,length(whichPts),length(loc_names));
bin_loc_nums = zeros(n_metrics,n_surrounds,length(whichPts),2);
output_for_plot = nan(length(whichPts),4,2);

for im = 1:n_metrics
    
    metric = all_metrics{im};
    
    for is = 1:n_surrounds
        
        surround = all_surrounds(is);
        
        %% prep table
        all_locs = array2table(nan(length(whichPts),length(loc_names)), 'VariableNames',loc_names);
        mt_vs_other = array2table(nan(length(whichPts),2), 'VariableNames',{'MT','all other except unspec'});
        
        for i = 1:length(whichPts)


            %% unchanegd labels, cblock
            cblock = out(i).change_block;
            unchanged_labels = out(i).unchanged_labels;
            switch metric
                case 'spike rate'
                    resp = out(i).rate;
                case 'node strength'
                    if do_norm
                        resp = out(i).metrics.ns_norm;
                    else
                        resp = out(i).metrics.ns;
                    end
            end
            
            %% Get anatomy
            unchanged_anatomy = out(i).unchanged_anatomy;

            %% Remove EKG and scalp electrodes
            ekg = identify_ekg_scalp(out(i).unchanged_labels);
            unchanged_labels(ekg) = [];
            resp(ekg,:) = [];
            unchanged_anatomy(ekg) = [];
            
            %% Define pre and post
            [pre,post] = get_surround_times(out(i).rate,cblock,surround);

     
            %% Group by localization
            [~,ana_loc] = anatomy_grouper(unchanged_anatomy);
            
            if 0
                % show result of anatomy grouper
                fprintf('\n%s:\n',out(i).name)
                table(unchanged_labels,unchanged_anatomy,ana_loc)
                pause
                
            end

            % loop through anatomy groups
            for j = 1:length(loc_names)

                % get the channels corresponding to that anatomy group
                chs = (strcmp(ana_loc,loc_names{j}));
                
                % get number of channels corresponding to that group
                loc_nums(im,is,i,j) = sum(chs);

                % Get the average rate change for those chs
                %% Get rate/ns relative change
                if strcmp(which_resp,'abs')
                    resp_change = (nanmean(resp(chs,post),'all') - nanmean(resp(chs,pre),'all'));
                elseif strcmp(which_resp,'rel')
                   % resp_change = (nansum(resp(chs,post),'all') - nansum(resp(chs,pre),'all'))./abs(nansum(resp(chs,pre),'all'));
                   % nansum slightly different. I think nanmean makes more
                   % sense because then, if one electrode disconnected, the
                   % nansum would essentially say it had fewer spikes
                   % because it does not divide by the number of non-nan
                   % channels. But results basically the same
                    resp_change = (nanmean(resp(chs,post),'all') - nanmean(resp(chs,pre),'all'))./abs(nanmean(resp(chs,pre),'all'));
                end
                avg_change = resp_change;
                
                %{
                if strcmp(which_resp,'abs')
                    resp_change = (nanmean(resp(chs,post),2) - nanmean(resp(chs,pre),2));
                elseif strcmp(which_resp,'rel')
                    resp_change = (nanmean(resp(chs,post),2) - nanmean(resp(chs,pre),2))./abs(nanmean(resp(chs,pre),2));
                end
                avg_change = nanmean(resp_change(chs));
                %}
                
                % add to locs table
                all_locs(i,j) = num2cell(avg_change);
                
                
                if 0
                    turn_nans_white(resp)
                    hold on
                    yticks(find(chs));
                    plot([cblock cblock],ylim,'r--')
                    xticks([pre(1) pre(end) post(1) post(end)])
                    colorbar
                    %plot([pre(1) pre(end)],[0 0],'g','linewidth',2)
                    %plot([post(1) post(end)],[0 0],'r','linewidth',2)
                    pause
                    
                end
            end
            
            %% the mt vs other analysis
            % Get those in MT group
            chs = strcmp(ana_loc,'mesial temporal');
            bin_loc_nums(im,is,i,1) = sum(chs);
            %avg_change = nanmean(resp_change(chs));
            %mt_vs_other(i,1) = num2cell(avg_change);
            
            
            % Get those in other groups (ignoring unspecified)
            chs = strcmp(ana_loc,'white matter') | ...
                strcmp(ana_loc,'temporal neocortical') | ... 
                strcmp(ana_loc,'other');
            bin_loc_nums(im,is,i,2) = sum(chs);
            %avg_change = nanmean(resp_change(chs));
            %mt_vs_other(i,2) = num2cell(avg_change);
            

        end
        
        % Convert to array
        all_locs_array = table2array(all_locs);
        
        % for Friedman, remove unspecified category
        all_locs_minus_unspecified = all_locs_array(:,[2:end]);
        
        % revmoe missing pt (none missing anymore so this does nothing)
        %{
        missing = nansum(all_locs_minus_unspecified,2) == 0;
        all_locs_minus_unspecified(missing,:) = [];
        %}
        
        if strcmp(which_test,'skillmack')
            [p stats] = skillmack(all_locs_minus_unspecified,1);
        elseif strcmp(which_test,'friedman')
            [p,~,stats] = friedman(all_locs_minus_unspecified,1,'off');
        end
        %
        if sum(isnan(all_locs_minus_unspecified)) == 0
            pf = friedman(all_locs_minus_unspecified,1,'off');
            if abs(p-pf) > 0.01
                error('what');
            end
        end
        %}
        chi2 = stats.T;
        df = stats.df;
        
        % make a string out of this
        curr_string = sprintf('chi^2(%d) = %1.1f, %s',df,chi2,pretty_p_text(p));
        final_array{im,is} = (curr_string);
        
        % Do between group analysis if p < 0.05/nsurrounds (this doesn't
        % happen)
        if p < 0.05/length(all_surrounds)
            
            if 0
                figure
                for k = 1:size(all_locs_minus_unspecified,2)
                    plot(k+rand(size(all_locs_minus_unspecified,1),1)*0.05,all_locs_minus_unspecified(:,k),'o')
                    hold on

                end
                xticks(1:size(all_locs_minus_unspecified,2))
                xticklabels(loc_names(2:end))
                title(all_metrics{im})
            end
            
            alpha_adjust = nchoosek(4,2); %6
            fprintf('\nFor a surround duration of %d:\n',all_surrounds(is));

            
            % loop over non-unspecified anatomy locs
            for j = 1:4
                
                % loop over the other anatomy locs
                for k = 1:j-1
                    
                    % perform a sign rank test (non parametric paired test)
                    % with alpha = 0.05/alpha_adjust (which is 6) to
                    % compare the two subgroups
                    alpha = 0.05/alpha_adjust;
                    psub = signrank(all_locs_minus_unspecified(:,j),all_locs_minus_unspecified(:,k));
                    if psub < alpha
                        fprintf(['\nPost-hoc analyses demonstrated a significant difference in relative %s change between '...
                            '%s (M (SD) %1.1f (%1.1f)) and %s (M (SD) %1.1f (%1.1f)) '...
                            '(Wilcoxon sign-rank, %s).\n'],all_metrics{im},...
                            loc_names{j+1},nanmean(all_locs_minus_unspecified(:,j)),...
                            nanstd(all_locs_minus_unspecified(:,j)),...
                            loc_names{k+1},nanmean(all_locs_minus_unspecified(:,k)),...
                            nanstd(all_locs_minus_unspecified(:,k)),pretty_p_text(psub));
                    else
                        fprintf(['\nThere was no significant difference between '...
                            '%s and %s'],loc_names{j+1},loc_names{k+1});
                    end
                    
                end
                
            end
        end
        
        % Show full result if main surround
        if is == main_surround
            if im == 1
                output_for_plot(:,:,1) = all_locs_minus_unspecified;
            elseif im == 2
               output_for_plot(:,:,2) = all_locs_minus_unspecified;
            end
            if strcmp(which_resp,'rel')
                fprintf(['\nThe mean (SD) relative %s change in the %d-hour peri-implant surround period was %1.1f (%1.1f)'...
                    ' for %s, %1.1f (%1.1f) for %s, %1.1f (%1.1f) for %s, and %1.1f (%1.1f) for %s regions.'...
                    ' The difference between groups was not significant (Skillings-Mack %s)\n'],...
                    all_metrics{im},all_surrounds(is),...
                    nanmean(all_locs_minus_unspecified(:,1)),nanstd(all_locs_minus_unspecified(:,1)),loc_names{2},...
                    nanmean(all_locs_minus_unspecified(:,2)),nanstd(all_locs_minus_unspecified(:,2)),loc_names{3},...
                    nanmean(all_locs_minus_unspecified(:,3)),nanstd(all_locs_minus_unspecified(:,3)),loc_names{4},...
                    nanmean(all_locs_minus_unspecified(:,4)),nanstd(all_locs_minus_unspecified(:,4)),loc_names{5},...
                    curr_string);
            elseif strcmp(which_resp,'abs')
                fprintf(['\nThe mean (SD) %s change in the %d-hour peri-implant surround period was %1.1f (%1.1f)'...
                    ' for %s, %1.1f (%1.1f) for %s, %1.1f (%1.1f) for %s, and %1.1f (%1.1f) for %s regions.'...
                    ' The difference between groups was not significant (Skillings-Mack %s)\n'],...
                    all_metrics{im},all_surrounds(is),...
                    nanmean(all_locs_minus_unspecified(:,1)),nanstd(all_locs_minus_unspecified(:,1)),loc_names{2},...
                    nanmean(all_locs_minus_unspecified(:,2)),nanstd(all_locs_minus_unspecified(:,2)),loc_names{3},...
                    nanmean(all_locs_minus_unspecified(:,3)),nanstd(all_locs_minus_unspecified(:,3)),loc_names{4},...
                    nanmean(all_locs_minus_unspecified(:,4)),nanstd(all_locs_minus_unspecified(:,4)),loc_names{5},...
                    curr_string);
                
            end
            
            % Plot it for my own understanding
            if 0
            figure
            for k = 1:size(all_locs_minus_unspecified,2)
                plot(k+rand(size(all_locs_minus_unspecified,1),1)*0.05,all_locs_minus_unspecified(:,k),'o')
                hold on

            end
            title(all_metrics{im})
            end
        end
        
        %% Do analysis for MT vs other
        %{
        mt_vs_other = table2array(mt_vs_other);
        [~,p,~,stats] = ttest(mt_vs_other(:,1),mt_vs_other(:,2));
        if is == main_surround
            fprintf(['\nThe mean (SD) relative %s change in the %d-hour peri-implant '...
                'period was %1.1f (%1.1f) for %s and %1.1f (%1.1f) for %s '...
                '(paired t-test: t(%d) = %1.1f, %s).\n'],...
                all_metrics{im},all_surrounds(is),...    
                nanmean(mt_vs_other(:,1)),nanstd(mt_vs_other(:,1)),loc_names{2},...
                nanmean(mt_vs_other(:,1)),nanstd(mt_vs_other(:,1)),'other localizations',...
                stats.df,stats.tstat,pretty_p_text(p));
            
            % Plot it for my own understanding
            if 0
            figure
            for k = 1:2
                plot(k+rand(size(mt_vs_other,1),1)*0.05,mt_vs_other(:,k),'o')
                hold on

            end
            title(all_metrics{im})
            end
        end
        %}
        
    end
end

%% show numbers
out_loc_names = loc_names(2:5);
s = main_surround;
im = 1;
    
curr_loc_nums = squeeze(loc_nums(im,s,:,:));

% Remove missing patient
curr_loc_nums(sum(curr_loc_nums(:,2:end),2)==0,:) = [];
%curr_loc_nums(:,1) = [];

% sum across patients
summed_loc_nums = sum(curr_loc_nums,1);
total_num = sum(summed_loc_nums);

% Count number of electrodes in each category
fprintf(['\n\nAcross all patients, %d (%1.1f%%) of electrode contacts were in %s,'...
    ' %d (%1.1f%%) were in %s regions, %d (%1.1f%%) were in %s regions,'...
    ' and %d (%1.1f%%) were in %s locations. %d (%1.1f%%) had an unspecified location.\n'],...
    summed_loc_nums(2),summed_loc_nums(2)/total_num*100,loc_names{2},...
    summed_loc_nums(3),summed_loc_nums(3)/total_num*100,loc_names{3},...
    summed_loc_nums(4),summed_loc_nums(4)/total_num*100,loc_names{4},...
    summed_loc_nums(5),summed_loc_nums(5)/total_num*100,loc_names{5},...
    summed_loc_nums(1),summed_loc_nums(1)/total_num*100);

numT = array2table(curr_loc_nums,'VariableNames',loc_names(1:5));
    


final_table = cell2table(final_array',...
    'RowNames',arrayfun(@(x) sprintf('%d',x),all_surrounds,...
    'UniformOutput',false),'VariableNames',all_metrics);


%writetable(final_table,[main_spike_results,'anatomy.csv'],'WriteRowNames',true)  

%% Also add to main supplemental table
writetable(final_table,[main_spike_results,'Supplemental Table 3.xlsx'],'Range','G2:H12','WriteVariableNames',false)


end