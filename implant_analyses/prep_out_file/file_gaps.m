function buffer = file_gaps(name)

switch name
    
    case 'HUP099'
        gap = '24:15:01';
    case 'HUP100'
        gap = '32:47:23';
    case 'HUP111'
        gap = '30:43:06';
    case 'HUP128'
        gap = '20:47:40';
    case 'HUP132'
        gap = '1:05:57';
    case 'HUP136'
        gap = '0:00:55';
    case 'HUP152'
        gap = '23:48:13';
    case 'HUP193'
        gap = '18:55:46';
    case 'HUP201'
        gap = '0:00:29';
    case 'HUP209'
        gap = '24:00:36';
    
end

% convert to duration
d = duration(str2double(strsplit(gap, ':')));

% convert to hours
h = hours(d);

% convert to blocks
buffer = round(h); 
% I am multiplying by 2 to get from hours to 30 minute blocks. But then I
% divide by 2 because I will add this buffer to both pre and post


end