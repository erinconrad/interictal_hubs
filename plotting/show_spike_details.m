function show_spike_details(values,HFdata,chLabels,details,fs,dur,filter)

%% basic info
gdf = details.filter(filter).gdf;
nsp = size(details.filter(filter).gdf,1);
sp_ch = gdf(:,1);
peak = details.filter(filter).peak;
rise = details.filter(filter).mid_rise;
fall = details.filter(filter).mid_fall;
unique_sp_chs = unique(sp_ch);
nchs = size(values,2);
n_un_chs = length(unique_sp_chs);

%% initialize figure and plotting things
figure
set(gcf,'position',[62 104 1145 701])

%% Plot the spikes, true data
%subplot(2,1,1)
offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for i = 1:n_un_chs 
    ich = unique_sp_chs(i);
    plot(values(:,ich)-offset,'k');
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
    plot(index,values(round(index),ch) - offset_sp,'bo')
    
    if ~isnan(rindex) && ~isnan(findex)
        plot(rindex,values(round(rindex),ch) - offset_sp,'go')
        plot(findex,values(round(findex),ch) - offset_sp,'ro')
    end
end


%% Plot HFdata
%{
subplot(2,1,2)
offset = 0;
ch_offsets = zeros(nchs,1);
ch_bl = zeros(nchs,1);

for i = 1:n_un_chs 
    
    ich = unique_sp_chs(i);
    
    
    
    plot(HFdata(:,ich)-offset,'k');
    ch_offsets(ich) = offset;
    ch_bl(ich) = -offset + nanmedian(HFdata(:,ich));
    hold on
    text(dur+1,ch_bl(ich),sprintf('%s',chLabels{ich}))
    
    if i<n_un_chs 
        if ~isnan(min(HFdata(:,ich)) - max(HFdata(:,unique_sp_chs(i+1))))
            offset = offset - (min(HFdata(:,ich)) - max(HFdata(:,unique_sp_chs(i+1))));
        end
    end
end

for s = 1:size(peak,1)
    index = peak(s);
    time = index/fs;
    
    rindex = rise(s);
    rtime = rise(s)/fs;
    
    findex = fall(s);
    ftime = fall(s)/fs;
    
    ch = sp_ch(s);
    offset_sp = ch_offsets(ch);

    value_sp = HFdata(round(index),ch);
    plot(index,value_sp - offset_sp,'bo')
    if ~isnan(rindex) && ~isnan(findex)
        plot(rindex,HFdata(round(rindex),ch) - offset_sp,'go')
        plot(findex,HFdata(round(findex),ch) - offset_sp,'ro')
    end
    
end
%}

pause
close(gcf)



end