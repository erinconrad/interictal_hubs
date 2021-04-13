function bad = reject_bad_chs(values,which_chs,chLabels,fs)

bad = [];

%% Parameters to reject super high variance
tile = 95;
mult = 10;
num_above = 5;

%% Parameter to reject high 60 Hz
percent_60_hz = 0.5;

for i = 1:length(which_chs)
    
    bad_ch = 0;
    
    ich = which_chs(i);
    eeg = values(:,ich);
    bl = median(eeg);
    
    %% Remove channels with any nans
    if sum(isnan(eeg)) > 0
        bad = [bad;ich];
        continue;
    end
    
    
    
    %% Remove channels if there are rare cases of super high variance above baseline (disconnection, moving, popping)
    pct = prctile(eeg,[100-tile tile]);
    thresh = [bl - mult*(bl-pct(1)), bl + mult*(pct(2)-bl)];
    sum_outside = sum(eeg > thresh(2) | eeg < thresh(1));
    if sum_outside > num_above
        bad_ch = 1;
    end
    
    if 0
        if bad_ch == 1
            plot(eeg,'r');
        else
            plot(eeg,'k');
        end
        title(chLabels{ich})
        hold on
        plot(xlim,[bl bl])
        plot(xlim,[thresh(1) thresh(1)]);
        plot(xlim,[thresh(2) thresh(2)]);
        title(sprintf('Sum outside: %d',sum_outside));
    end
    
    if bad_ch == 1
        bad = [bad;ich];
        continue;
    end
    
    %% Remove channels with a lot of 60 Hz noise, suggesting poor impedance
    
    % Calculate fft
    Y = fft(eeg-mean(eeg));
    
    % Get power
    P = abs(Y).^2;
    freqs = linspace(0,fs,length(P)+1);
    freqs = freqs(1:end-1);
    
    % Take first half
    P = P(1:ceil(length(P)/2));
    freqs = freqs(1:ceil(length(freqs)/2));
    
    P_60Hz = sum(P(freqs > 58 & freqs < 62))/sum(P);
    if P_60Hz > percent_60_hz
        bad_ch = 1;
    end
    
    if 0
        spectrogram(eeg,[],[],[],fs);
        title(sprintf('Percent 60 Hz power %1.1f',P_60Hz*100))
    end
    
    if bad_ch == 1
        bad = [bad;ich];
        continue;
    end
    
end

end