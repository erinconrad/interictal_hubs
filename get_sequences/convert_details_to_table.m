function T = convert_details_to_table(spikes,which_files,which_blocks)

orig_which_blocks = which_blocks;

if isempty(which_files)
    which_files = 1:length(spikes.file);
end


%% Initialize stuff
ch = [];
peak_idx = [];
peak_amp = [];
first_rise_idx = [];
mid_rise_idx = [];
mid_fall_idx = [];
last_fall_idx = [];
fwhm = [];
wfile = [];
run_start = [];
fs = [];


for i = 1:length(which_files)
    f = which_files(i);
    
    file = spikes.file(f);
    
    if isempty(orig_which_blocks)
        which_blocks = 1:length(file.block);
    end
    
    for j = 1:length(which_blocks)
        
        
        
        b = which_blocks(j);
        block = file.block(b);
        
        if isempty(block.details)
            continue;
        end
        ns = length(block.details.ch);
        ch = [ch;block.details.ch];
        peak_idx = [peak_idx;block.details.peak_idx];
        peak_amp = [peak_amp;block.details.peak_amp];
        first_rise_idx = [first_rise_idx;block.details.first_rise_idx];
        mid_rise_idx = [mid_rise_idx;block.details.mid_rise_idx]; 
        mid_fall_idx = [mid_fall_idx;block.details.mid_fall_idx];
        last_fall_idx = [last_fall_idx;block.details.last_fall_idx];
        fwhm = [fwhm;block.details.fwhm_ins];
        wfile = [wfile;repmat(f,ns,1)];
        run_start = [run_start;repmat(block.run_times(1),ns,1)];
        fs = [fs;repmat(block.fs,ns,1)];
        
    end
    
end



T = table(wfile,run_start,fs,ch,peak_idx,peak_amp,first_rise_idx,mid_rise_idx,...
    mid_fall_idx,last_fall_idx,fwhm);


end