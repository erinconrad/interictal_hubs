function details = get_spike_details(gdf,values,fs)

%% Parameters
time_to_look = [-50e-3 50e-3];
perc_rise = 0.1;
look_before_peak = -200e-3;
tmul = 3;
min_rise_time = 10e-3;
fn_fr  = 7; % high pass filter for spikey component
fr     = 40; % low pass filter for spikey component
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
    orig_s_index = gdf(s,3);
    
    data = values(:,ch);
    fndata   = eegfilt(data, fn_fr, 'hp',fs); % high pass filter
    HFdata    = eegfilt(fndata, fr, 'lp',fs); % low pass filter
    
    %% Find spike peak time and amplitude
    % Restrict times to look 
    idx_to_look = max(1,round(t+time_to_look(1)*fs)):min(round(t+time_to_look(2)*fs),length(HFdata));

    % Find the peak in the data
    [peak_amp,ind] = max(abs(HFdata(idx_to_look)));
    

    % re-align index
    abs_ind = ind + idx_to_look(1) - 1;

    if 0
        plot(data)
        hold on
        plot(HFdata)
        plot(abs_ind,data(abs_ind),'o')
    end

    peak_dev = HFdata(abs_ind); % this is signed
    thresh_dev = tmul * std(HFdata);
    
    %% Fill peak amp and peak time
    details.peak_amp(s) = abs(HFdata(abs_ind));
    details.peak_idx(s) = abs_ind;
    details.ch(s) = ch;
    
    %% Find initial rise and mid-rise time

    % Look before the peak and find the first point at which all subsequent
    % points leading up to the peak are a minimum deviation from baseline
 
    idx_to_peak = max(1,round(abs_ind+look_before_peak*fs)):abs_ind;
    rise = abs(HFdata(idx_to_peak)) > thresh_dev;
    num_rise = cumsum(rise);
    min_rise_idx = min_rise_time * fs;
    all_enough_rise = find(num_rise > min_rise_idx);
    
    if isempty(all_enough_rise)
        
        details.first_rise_idx(s) = nan;
        details.last_fall_idx(s) = nan;
        details.mid_rise_idx(s) = nan;
        details.mid_fall_idx(s) = nan;
        details_fwhm_ins(s) =  nan;
        continue
    end
    
    first_rise = all_enough_rise(1) - min_rise_idx + idx_to_peak(1) - 1;

    if 0 
    plot(HFdata)
    hold on
    plot(idx_to_peak(rise),HFdata(idx_to_peak(rise)))
    plot(xlim,[thresh_dev thresh_dev])
    plot(xlim,-[thresh_dev thresh_dev])
    plot(abs_ind,HFdata(abs_ind),'o')
    plot(first_rise,HFdata(first_rise),'x')
    end

    mid_rise_time = round(mean([abs_ind,first_rise]));

    %% Find fall
    idx_after_peak = abs_ind:min(round(abs_ind-look_before_peak*fs),length(HFdata));
    rise = abs(HFdata(idx_after_peak)) > thresh_dev;
    num_rise = cumsum(rise,'reverse');
    min_rise_idx = round(min_rise_time * fs);
    all_enough_rise = find(num_rise > min_rise_idx);
    
    if isempty(all_enough_rise)
        
        details.first_rise_idx(s) = nan;
        details.last_fall_idx(s) = nan;
        details.mid_rise_idx(s) = nan;
        details.mid_fall_idx(s) = nan;
        details_fwhm_ins(s) =  nan;
        continue
        
        
    end
    
    last_rise = idx_after_peak(end) - (length(idx_after_peak)-all_enough_rise(end)) + min_rise_idx + 1;
    mid_fall_time = round(mean([abs_ind,last_rise]));


    if 0
    figure
    plot(HFdata)
    hold on
    plot(idx_after_peak(rise),HFdata(idx_after_peak(rise)))
    plot(xlim,[thresh_dev thresh_dev])
    plot(xlim,-[thresh_dev thresh_dev])
    plot(abs_ind,HFdata(abs_ind),'o')
    plot(first_rise,HFdata(first_rise),'x')
    plot(last_rise,HFdata(last_rise),'+')
    end


    if 0
        figure
        plot(idx_to_peak(1):idx_after_peak(end),data(idx_to_peak(1):idx_after_peak(end)));
        hold on
        plot(abs_ind,data(abs_ind),'bo')
        plot(mid_rise_time,data(mid_rise_time),'g+')
        plot(mid_fall_time,data(mid_fall_time),'r+')
        pause
        close(gcf)
    end


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