function details = new_spike_details(gdf,values,hf_values,fs)

%% Parameters
peak_look = [-50e-3 50e-3];
min_look = [-100e-3 100e-3];
ns = size(gdf,1);
fn_fr  = 7;
fr     = 40;

%% Initialize structures
%details.orig_gdf = gdf;
for f = 1:2
    details.filter(f).peak = nan(ns,1);
    details.filter(f).amp = nan(ns,1);
    details.filter(f).rise = nan(ns,1);
    details.filter(f).fall = nan(ns,1);
    details.filter(f).gdf = gdf;
end

% Loop over spikes
for s = 1:size(gdf,1)
    
    ch = gdf(s,1);
    t = gdf(s,2);
    
    orig_data = values(:,ch)-nanmean(values(:,ch));
    HF_data    = hf_values(:,ch); 
    %HF_data = eegfilt(orig_data, fn_fr, 'hp',fs);
    %HF_data    = eegfilt(HF_data, fr, 'lp',fs);
    
    % Loop over methods of filtering (HF filter or no)
    for f = 1:2
        if f == 1
            data = orig_data;
        else
            data = HF_data;
        end
        
        %% Find peak
        % Restrict times to look 
        idx_to_look = max(1,round(t+peak_look(1)*fs)):min(round(t+peak_look(2)*fs),length(data));

        % Find the peak in the data
        [peak_amp,ind] = max(abs(data(idx_to_look))); % unsigned
        
        
        % re-align index
        abs_ind = ind + idx_to_look(1) - 1;
        
        % Add to struct
        details.filter(f).peak(s) = abs_ind;
        details.filter(f).amp(s) = peak_amp;
        
        %% Find initial rise and mid-rise time
    
        
        % Define times to look (from 100 ms before peak to peak)
        idx_to_peak = max(1,round(abs_ind+min_look(1)*fs)):abs_ind;
        
        % if peak dev is positive, look for min, else look for max
        signed_peak = data(abs_ind);
        if signed_peak > 0
            [~,min_before_peak] = min(data(idx_to_peak));
        else
            [~,min_before_peak] = max(data(idx_to_peak));
        end
        
        % time of this initial rise
        first_rise = min_before_peak + idx_to_peak(1) - 1;

        % mid rise time (halfway from first rise to peak)
        mid_rise_time = round(mean([abs_ind,first_rise]));
        
        % Add to struct
        details.filter(f).rise(s) = first_rise;
        
        %% Find fall
        idx_after_peak = abs_ind:min(round(abs_ind+min_look(2)*fs),length(data));
        % if peak dev, positive, look for min, else look for max
        if signed_peak > 0
            [~,min_after_peak] = min(data(idx_after_peak));
        else
            [~,min_after_peak] = max(data(idx_after_peak));
        end

        last_rise = min_after_peak + idx_after_peak(1) - 1;
        mid_fall_time = round(mean([abs_ind,last_rise]));
        
        % Add to struct
        details.filter(f).fall(s) = last_rise;
        
        %% Plot
        if 0
            figure
            plot(data)
            hold on
            plot(abs_ind,data(abs_ind),'bo')
            hold on
            plot(first_rise,data(first_rise),'go')
            plot(last_rise,data(last_rise),'ro')
            plot(mid_rise_time,data(mid_rise_time),'g+')
            plot(mid_fall_time,data(mid_fall_time),'r+')
            pause
            close(gcf)
            
        end
        
    end
    
    
    
end


%% Now loop over pairs of spikes in the gdf and remove whichever one has the lowest amplitude
% Loop over gdf, 2 at a time
%{
keep = ones(ns,1);
for s = 1:2:ns-1
    amp1 = details.filter(1).amp(s);
    amp2 = details.filter(1).amp(s+1);

    if amp1 >= amp2
        keep(s+1) = 0;
    else
        keep(s) = 0;
    end
end


for f = 1:2
    details.filter(f).amp(keep==0) = [];
    details.filter(f).peak(keep==0) = [];
    details.filter(f).rise(keep==0) = [];
    details.filter(f).fall(keep==0) = [];
    details.filter(f).gdf(keep==0,:) = [];
    details.filter(f).gdf(:,3) = [];
end
%}

%{
for f = 1:2
    % Loop over gdf, 2 at a time
    keep = ones(ns,1);
    for s = 1:2:ns-1
        amp1 = details.filter(f).amp(s);
        amp2 = details.filter(f).amp(s+1);

        if amp1 >= amp2
            keep(s+1) = 0;
        else
            keep(s) = 0;
        end
    end

    details.filter(f).amp(keep==0) = [];
    details.filter(f).peak(keep==0) = [];
    details.filter(f).rise(keep==0) = [];
    details.filter(f).fall(keep==0) = [];
    details.filter(f).gdf(keep==0,:) = [];
    details.filter(f).gdf(:,3) = [];

%}
end
