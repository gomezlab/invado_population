function find_cell_mask_properties_experiment(base_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('gelatin_min_value',382,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

addpath('matlab_scripts');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

for i=1:length(fields)
    try
        find_cell_mask_properties(fullfile(base_dir,fields(i).name), ...
            'gelatin_min_value',i_p.Results.gelatin_min_value);
        disp(['Done with ', fullfile(base_dir,fields(i).name)]);
    catch %#ok<CTCH>
        disp(['Problem with ', fullfile(base_dir,fields(i).name)]);
    end
end

toc;