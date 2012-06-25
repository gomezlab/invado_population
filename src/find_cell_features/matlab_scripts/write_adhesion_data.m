function write_adhesion_data(S,varargin)
% WRITE_STRUCT_DATA     write most the data stored in a given struct to a
%                       set of ascii formated files, using the field
%                       names as the file names
%
%   write_struct_data(S) writes most the field names in struct 'S' to ascii
%   formated files, suitable for use in other programs, the fieldnames are
%   used as the file names, all of the files are placed in subfolder
%   'raw_data' of the current working directory
%
%   write_struct_data(S,'out_dir',d) writes most the field names in struct
%   'S' to ascii formated files, suitable for use in other programs, the
%   fieldnames are used as the file names, all of the files are placed in
%   the provided directory 'd'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'WRITE_ADHESION_DATA';

i_p.addRequired('S',@isstruct);
i_p.addParamValue('out_dir','raw_data',@ischar);

i_p.parse(S,varargin{:});
out_dir = i_p.Results.out_dir;

if (not(exist(out_dir,'dir')))
    mkdir(out_dir);
end

to_exclude = {'ConvexHull','ConvexImage','Image','FilledImage', ...
    'PixelList', 'SubarrayIdx', 'Border_pix', 'Extrema'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

field_names = fieldnames(S);

ad_props_cell = struct2cell(S);

print_strings = struct('PixelIdxList','%0.f','Cell_gel_diff','%0.10f', ... 
    'Cell_gel_diff_p_val','%0.10f','Overlap_region_size','%d');

for i = 1:size(field_names,1)
    
    if(strmatch(field_names(i),to_exclude))
        continue;
    end
    
    format_string = '%f';
    if(strmatch(field_names(i),fieldnames(print_strings)))
        format_string = print_strings.(field_names{i});
    end    
    
    file_out = fullfile(out_dir,[cell2mat(field_names(i)),'.csv']);
    
    data = ad_props_cell(i,:);
    if (any(size(data) == 0))
        continue;
    end
    output_CSV_from_cell(data, file_out, 'format', format_string);
end

