function dir_ls = filter_to_time_series(dir_ls)
%FILTER_TO_TIME_SERIES    Takes in a dir struct, outputs the same dir
%                         struct with all the folder names that includes
%                         'time_series' 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

records_to_keep = ones(length(dir_ls),1);
for i=1:length(dir_ls)
    if (isempty(findstr(dir_ls(i).name,'time_series')))
        records_to_keep(i) = 0;
    end
end

dir_ls = dir_ls(logical(records_to_keep));