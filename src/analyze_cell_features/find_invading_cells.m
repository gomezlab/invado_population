function find_invading_cells(field_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('field_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('median_filter',-1.5,@isnumeric);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(field_dir,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_data = struct();
files = struct();
files.tracking = fullfile(field_dir,'tracking_matrices','tracking_seq.csv');
%check for the tracking file first, if absent, no cells were found, return
%from the function with a message, but no error
if (not(exist(files.tracking,'file')))
    disp('No tracking matrix found, assuming no cells in field');
    return;
end

data_series_folder = fullfile(field_dir,'cell_props','lin_time_series');

files.p_vals = fullfile(data_series_folder,'Cell_gel_diff_p_val.csv');
files.cell_diffs = fullfile(data_series_folder,'Cell_gel_diff.csv');
files.cell_diff_medians = fullfile(data_series_folder,'Cell_gel_diff_median.csv');
files.cell_diff_percents = fullfile(data_series_folder,'Cell_gel_diff_percent.csv');
files.gel_diff_minus = fullfile(data_series_folder,'Gel_diff_minus_surrounding.csv');

files.cell_total_degrade = fullfile(data_series_folder,'..','final_gel_diffs.csv');
files.corrected_cell_total_degrade = fullfile(data_series_folder,'..','corrected_final_gel_diffs.csv');

these_types = fieldnames(files);
for j = 1:length(these_types)
    this_file = files.(these_types{j});
    
    %matlab doesn't like you to reference fields that haven't been
    %created, so create files that aren't present yet before loading
    %data in
    if(isempty(strcmp(these_types{j},fieldnames(raw_data))))
        raw_data.(these_types{j}) = [];
    end
    
    if (exist(this_file,'file'))
        raw_data.(these_types{j}) = load(this_file);
    else
        error('Invado:MissingFile',['Can''t find ',this_file])
    end
end

raw_data.cell_total_degrade = repmat(raw_data.cell_total_degrade,1,size(raw_data.cell_diff_percents,2));
raw_data.corrected_cell_total_degrade = repmat(raw_data.corrected_cell_total_degrade,1,size(raw_data.cell_diff_percents,2));

%check that all the raw data files are the same size
these_names = fieldnames(raw_data);
poss_name_combinations = combnk(1:length(these_names),2);
for j=1:size(poss_name_combinations,1)
    assert(all(size(raw_data.(these_names{poss_name_combinations(j,1)})) == ...
        size(raw_data.(these_names{poss_name_combinations(j,2)}))))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process and Output 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
processed_data = process_raw_data(raw_data);

output_dir = fullfile(field_dir,'cell_props');

csvwrite(fullfile(output_dir,'active_degrade.csv'),processed_data.active_degrade);

csvwrite(fullfile(output_dir,'longevity.csv'),processed_data.longevities);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function process_data = process_raw_data(raw_data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('raw_data',@isstruct);
i_p.addParamValue('filter_set',NaN,@islogical);

i_p.parse(raw_data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%use the filter_set variable to filter all the data sets before continuing
if (isempty(strcmp('filter_set',i_p.UsingDefaults)))
    these_names = fieldnames(raw_data);
    for j=1:size(these_names,1)
        raw_data.(these_names{j}) = raw_data.(these_names{j})(i_p.Results.filter_set,:);
    end
end

process_data = struct();

process_data.active_degrade = not(isnan(raw_data.gel_diff_minus)) & raw_data.gel_diff_minus < -1.6 ...
    & not(isnan(raw_data.corrected_cell_total_degrade)) & raw_data.corrected_cell_total_degrade < -0.4;

disp(['Detected ', num2str(sum(process_data.active_degrade(:))), ' invasion events.']);
% disp(['Bonferroni Corrected p-value threshold: ', num2str(0.05/bonferroni_correction)]);

process_data.live_cells = raw_data.tracking > -1;
process_data.longevities = sum(process_data.live_cells,2)/2;

process_data.ever_degrade = [];
for i=1:size(raw_data.tracking,1)
    process_data.ever_degrade = [process_data.ever_degrade, any(process_data.active_degrade(i,:))];
end

process_data.has_degraded = zeros(size(raw_data.tracking));
for i=1:size(raw_data.tracking,1)
    for j = 1:size(raw_data.tracking,2)
        process_data.has_degraded(i,j) = any(process_data.has_degraded(i,1:j));
    end
end

process_data.degrade_percentage = sum(process_data.has_degraded)/size(raw_data.tracking,1);