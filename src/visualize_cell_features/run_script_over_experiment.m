function run_script_over_experiment(base_dir,script,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('script',@(x)ischar(x));

i_p.addParamValue('gelatin_min_value',382,@(x)isnumeric(x) && x > 0);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,script,varargin{:});

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

for i=1:length(fields)
    try
        execution_string = [script,'(''',fullfile(base_dir,fields(i).name),''');'];
        eval(execution_string);
        disp(['Done with ', fullfile(base_dir,fields(i).name)]);
    catch %#ok<CTCH>
        disp(['Problem with ', fullfile(base_dir,fields(i).name)]);
    end
end

toc;