function values = pull_ieeg_data(fname, login_name, pwfile, run_idx)

session = IEEGSession(fname, login_name, pwfile);
channelLabels = session.data.channelLabels;
nchs = size(channelLabels,1);

% Break the number of channels in half to avoid wacky server errors
values1 = session.data.getvalues(run_idx,1:floor(nchs/4));
values2 = session.data.getvalues(run_idx,floor(nchs/4)+1:floor(2*nchs/4));
values3 = session.data.getvalues(run_idx,floor(2*nchs/4)+1:floor(3*nchs/4)); 
values4 = session.data.getvalues(run_idx,floor(3*nchs/4)+1:nchs); 

session.delete;

values = [values1,values2,values3,values4];

end