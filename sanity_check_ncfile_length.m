% sanity check to make sure the g16 sst is correct:
basetimenum = datenum('2020-01-01','yyyy-mm-dd');

days_of_year = 1:59;
data_root = '/Volumes/sst_raw/g16';


for id = 1:length(days_of_year)
    
    DOY = days_of_year(id);
    path2data = [data_root filesep num2str(DOY, '%3.3i') filesep 'l3c'];
    dataID_date = datestr(basetimenum + (DOY-1),'yyyymmdd');
    
      
    % for it = 1 :nhr
    ncfiles = dir([path2data filesep '*.nc']);
    dataFNs = {ncfiles.name};
    
    clear GOES_ATOMIC
    disp([dataID_date ', nt = ' num2str(length(dataFNs))]);
end