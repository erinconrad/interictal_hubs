function show_spike_details(values,chLabels,details,fs,dur)

%% Parameters
fn_fr  = 7; % high pass filter for spikey component
fr     = 40; % low pass filter for spikey component

%% basic info
gdf = details.gdf;
nsp = size(details.gdf,1);
sp_ch = details.ch;
peak = details.peak_idx;
rise = details.mid_rise_idx;
fall = details.mid_fall_idx;
unique_sp_chs = unique(gdf(:,1));
nchs = size(values,2);
n_un_chs = length(unique_sp_chs);

%% initialize figure and plotting things
figure
set(gcf,'position',[62 104 1145 701])

%% Plot the spikes, true data
subplot(2,1,1)
offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for i = 1:n_un_chs 
    ich = unique_sp_chs(i);
    plot(linspace(0,dur,size(values,1)),values(:,ich)-offset,'k');
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(values(:,ich));
    hold on
    text(dur+0.05,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
    if i<n_un_chs 
        if ~isnan(min(values(:,ich)) - max(values(:,unique_sp_chs(i+1))))
            offset = offset - (min(values(:,ich)) - max(values(:,unique_sp_chs(i+1))));
        end
    end
end

for s = 1:size(peak,1)
    %index = spikes(s,1);
    index = peak(s);
    time = index/fs;
    
    rindex = rise(s);
    rtime = rise(s)/fs;
    
    findex = fall(s);
    ftime = fall(s)/fs;
    
    ch = sp_ch(s);
    offset_sp = ch_offsets(ch);
    plot(time,values(round(index),ch) - offset_sp,'bo')
    if ~isnan(rindex) && ~isnan(findex)
        plot(rtime,values(round(rindex),ch) - offset_sp,'go')
        plot(ftime,values(round(findex),ch) - offset_sp,'ro')
    end
end

%% Get HFdata 
HFdata = zeros(size(values));
for i = 1:size(values,1)
    data = values(:,ich);
    fndata   = eegfilt(data, fn_fr, 'hp',fs); % high pass filter
    HFdata(:,i) = eegfilt(fndata, fr, 'lp',fs); % low pass filter
    
end

%% Plot HFdata
subplot(2,1,2)
offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for i = 1:n_un_chs 
    
    ich = unique_sp_chs(i);
    
    
    
    plot(linspace(0,dur,size(values,1)),HFdata(:,ich)-offset,'k');
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(HFdata(:,ich));
    hold on
    text(dur+0.05,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
    if i<n_un_chs 
        if ~isnan(min(HFdata(:,ich)) - max(HFdata(:,unique_sp_chs(i+1))))
            offset = offset - (min(HFdata(:,ich)) - max(HFdata(:,unique_sp_chs(i+1))));
        end
    end
end

for s = 1:size(peak,1)
    %index = spikes(s,1);
    index = peak(s);
    time = index/fs;
    
    ch = sp_ch(s);
    offset_sp = ch_offsets(ch);

    value_sp = HFdata(round(index),ch);
    plot(time,value_sp - offset_sp,'ro')
    
end

pause
close(gcf)



end