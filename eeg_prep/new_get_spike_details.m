function details = new_get_spike_details(gdf,values,hf_values,fs)

%% Parameters
time_to_look = [-50e-3 50e-3];
look_before_peak = -100e-3;
tmul = 3;
ns = size(gdf,1);

%% Initialize structures of useful info
details.peak_amp = nan(ns,1);
details.peak_idx = nan(ns,1);
details.ch = nan(ns,1);
details.first_rise_idx = nan(ns,1);
details.last_fall_idx = nan(ns,1);
details.mid_rise_idx = nan(ns,1);
details.mid_fall_idx = nan(ns,1);
details.fwhm_ins = nan(ns,1);
details.gdf = gdf;


% Loop over spikes
for s = 1:size(gdf,1)
    
    ch = gdf(s,1);
    t = gdf(s,2);
    
    % subtract baseline
    data = values(:,ch)-median(values(:,ch));
    HFdata    = hf_values(:,ch); 
    
    %% Find spike peak time and amplitude
    % Restrict times to look 
    idx_to_look = max(1,round(t+time_to_look(1)*fs)):min(round(t+time_to_look(2)*fs),length(HFdata));

    % Find the peak in the data
    [~,ind] = max(abs(data(idx_to_look)));
    

    % re-align index
    abs_ind = ind + idx_to_look(1) - 1;

    if 0
        plot(data)
        hold on
        plot(abs_ind,data(abs_ind),'o')
        pause
        close(gcf)
    end

    peak_dev = data(abs_ind); % this is signed
    thresh_dev = tmul * std(HFdata);
    
    %% Fill peak amp and peak time
    details.peak_amp(s) = abs(data(abs_ind));
    details.peak_idx(s) = abs_ind;
    details.ch(s) = ch;
    
    %% Find initial rise and mid-rise time
    % Look back this far
    idx_to_peak = max(1,round(abs_ind+look_before_peak*fs)):abs_ind;
    
    % if peak dev is positive, look for min, else look for max
    if peak_dev > 0
        [~,min_before_peak] = min(data(idx_to_peak));
    else
        [~,min_before_peak] = max(data(idx_to_peak));
    end
    
    % time of this initial rise
    first_rise = min_before_peak + idx_to_peak(1) - 1;

    % mid rise time (halway from first rise to peak)
    mid_rise_time = round(mean([abs_ind,first_rise]));

    %% Find fall
    idx_after_peak = abs_ind:min(round(abs_ind-look_before_peak*fs),length(HFdata));
    % if peak dev, positive, look for min, else look for max
    if peak_dev > 0
        [~,min_after_peak] = min(data(idx_after_peak));
    else
        [~,min_after_peak] = max(data(idx_after_peak));
    end
    
    last_rise = min_after_peak + idx_after_peak(1) - 1;
    mid_fall_time = round(mean([abs_ind,last_rise]));


    if 0
        figure
        subplot(2,1,1)
        plot(data)
        hold on
        plot(abs_ind,data(abs_ind),'bo')
        hold on
        plot(first_rise,data(first_rise),'go')
        plot(last_rise,data(last_rise),'ro')
        plot(mid_rise_time,data(mid_rise_time),'g+')
        plot(mid_fall_time,data(mid_fall_time),'r+')

        subplot(2,1,2)
        plot(HFdata)
        hold on
        plot(abs_ind,HFdata(abs_ind),'bo')
        hold on
        plot(first_rise,HFdata(first_rise),'go')
        plot(last_rise,HFdata(last_rise),'ro')
        plot(mid_rise_time,HFdata(mid_rise_time),'g+')
        plot(mid_fall_time,HFdata(mid_fall_time),'r+')
        pause
        close(gcf)
    end

   
    details.first_rise_idx(s) = first_rise;
    details.last_fall_idx(s) = last_rise;
    details.mid_rise_idx(s) = mid_rise_time;
    details.mid_fall_idx(s) = mid_fall_time;
    details_fwhm_ins(s) = (mid_fall_time - mid_rise_time)/fs;
    
end

%% Now loop over pairs of spikes in the gdf and remove whichever one has the lowest amplitude

% Loop over gdf, 2 at a time
keep = ones(ns,1);
for s = 1:2:ns-1
    amp1 = details.peak_amp(s);
    amp2 = details.peak_amp(s);
    
    if amp1 >= amp2
        keep(s+1) = 0;
    else
        keep(s) = 0;
    end
end

details.peak_amp(keep==0) = [];
details.ch(keep==0) = [];
details.peak_idx(keep==0) = [];
details.first_rise_idx(keep==0) = [];
details.last_fall_idx(keep==0) = [];
details.mid_rise_idx(keep==0) = [];
details.mid_fall_idx(keep==0) = [];
details.fwhm_ins(keep==0) = [];
details.exp_gdf = details.gdf;
details.gdf(keep==0,:) = [];
details.gdf(:,3) = [];





end