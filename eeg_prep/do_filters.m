function out = do_filters(eeg,fs)

out = eeg;
for ich = 1:size(eeg,2)
    vals = eeg(:,ich);
  
    %% Turn nans into median of signal
    vals(isnan(vals)) = nanmedian(vals);
    
    %% Subtract the baseline
    vals = vals-nanmedian(vals);
    
    %% Notch filter
    vals = bandstop(vals,[59 61],fs,...
        'ImpulseResponse','iir');
    
    %% Bandpass filter
    d = designfilt('bandpassiir','FilterOrder',8, ...
    'HalfPowerFrequency1',1,'HalfPowerFrequency2',70, ...
    'SampleRate',fs);
    vals = filter(d,vals);

    
    out(:,ich) = vals;
    
    if 0
        figure
        plot(eeg(:,ich))
        hold on
        plot(vals)
        pause
        close(gcf)
        
    end
    
end

end