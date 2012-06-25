function output_CSV_from_cell(data, out_file, varargin)
% output_CSV_from_cell    writes a provided cell data structure to a CSV
%                         file
%
%   write_struct_data(D,OF) writes the data in cell 'D' to file 'DF', using
%   fprintf format '%f'
%
%   Options:
%       -format: specify the string passed to fprintf for number conversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;
i_p.FunctionName = 'output_CSV_from_cell';

i_p.addRequired('data',@iscell);
i_p.addRequired('out_file', @(x) exist(fileparts(x),'dir') == 7 );

i_p.addParamValue('format','%f',@ischar);

i_p.parse(data,out_file,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

file_handle = fopen(i_p.Results.out_file,'wt');
for i = 1:max(size(data))
    for j = 1:max(size(data{i}))
        assert(any(size(data{i}) == 1))
        
        if (size(data{i},1) == 0)
            continue;
        end
        
        if (j < max(size(data{i})))
            fprintf(file_handle,[i_p.Results.format,','],data{i}(j));
        else
            fprintf(file_handle,[i_p.Results.format,'\n'],data{i}(j));
        end
    end
end
fclose(file_handle);
