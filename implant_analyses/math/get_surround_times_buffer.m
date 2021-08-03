function [pre,post] = get_surround_times_buffer(rate,cblock,surround,buffer,...
    pre_nans,post_nans)

%{
This function takes the spike rate data, the time of the re-implant, and
returns times to call pre-implant and post-implant. Other inputs are the
"surround" time (how many hours pre and post to take) and the "buffer" (how
far before and after the reimplant time to ignore. The buffer is only used
in the Monte Carlo simulations. This is because for most patients, there is
a gap in ieeg recording around the implant revision time. To avoid a bias,
I need to include this gap in the Monte Carlo simulations.
%}


%% Establish the pre and post times, taking into account nans and buffer
pre = max(1,cblock-surround-buffer-pre_nans):...
    max(1,cblock-buffer-pre_nans); % go buffer back from the first non nan time
post = min(cblock+buffer+post_nans,size(rate,2)):...
    min(cblock+surround+buffer+post_nans,size(rate,2)); % go buffer ahead from the first non nan time


if 0
    figure
    turn_nans_white(rate)
    hold on
    plot([cblock cblock],ylim,'r--','linewidth',3)
    plot([pre(1) pre(1)],ylim,'g--','linewidth',3)
    plot([pre(end) pre(end)],ylim,'g--','linewidth',3)
    
    plot([post(1) post(1)],ylim,'b--','linewidth',3)
    plot([post(end) post(end)],ylim,'b--','linewidth',3)
    
    pause
    close(gcf)
end

end