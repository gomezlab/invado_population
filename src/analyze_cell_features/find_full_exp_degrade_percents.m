function find_full_exp_degrade_percents(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('min_longevity',10,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(exp_dir);
fields = filter_to_time_series(fields);

if (isempty(fields))
    disp('Expected to find directories for the fields, did you provide the correct directory?');
    return;
end

raw_data = struct('active_degrade',[],'longevity',[],'tracking',[], ...
    'gel_minus_surrounding',[],'corrected_final_diff',[],'field_number',[],...
    'cell_number',[],'cell_speed',[],'degradation_area',[]);

for j=1:length(fields)
    props_base = fullfile(exp_dir,fields(j).name,'cell_props');

    tracking_file = fullfile(exp_dir,fields(j).name,'tracking_matrices','tracking_seq.csv');
    if (not(exist(tracking_file,'file')))
        fprintf('No tracking file found in field %s, moving to next set.\n',fields(j).name);
        continue;
    end
    
    degrade_data = csvread(fullfile(props_base,'active_degrade.csv'));
    raw_data.active_degrade = [raw_data.active_degrade; degrade_data];
    
    longev_data = csvread(fullfile(props_base,'longevity.csv'));
    raw_data.longevity = [raw_data.longevity;longev_data];

    degradation_area = csvread(fullfile(props_base,'degradation_area.csv'));
    raw_data.degradation_area = [raw_data.degradation_area;degradation_area];
    
    gel_minus = csvread(fullfile(props_base,'lin_time_series','Gel_diff_minus_surrounding.csv'));
    raw_data.gel_minus_surrounding = [raw_data.gel_minus_surrounding;gel_minus];    

    cell_speed = csvread(fullfile(props_base,'lin_time_series','Cell_speed.csv'));
    raw_data.cell_speed = [raw_data.cell_speed;cell_speed];    
  
    corrected_final_diff = csvread(fullfile(props_base,'corrected_final_gel_diffs.csv'));
    raw_data.corrected_final_diff = [raw_data.corrected_final_diff;corrected_final_diff];    
    
    tracking_data = csvread(fullfile(exp_dir,fields(j).name,'tracking_matrices','tracking_seq.csv'));
    raw_data.tracking = [raw_data.tracking;tracking_data];
    
    raw_data.field_number = [raw_data.field_number;j*ones(size(tracking_data,1),1)];
    raw_data.cell_number = [raw_data.cell_number;(1:size(tracking_data,1))'];
end

longev_filter = raw_data.longevity > i_p.Results.min_longevity;

processed_data = process_raw_data(raw_data,'filter_set',longev_filter);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_dir = fullfile(exp_dir,'overall_results');
if (not(exist(output_dir,'dir')))
    mkdir(output_dir)
end

csvwrite(fullfile(output_dir,'degrade_percentage.csv'),processed_data.degrade_percentage);
csvwrite(fullfile(output_dir,'has_degraded.csv'),processed_data.has_degraded);

csvwrite(fullfile(output_dir,'fields_cells.csv'),...
     [processed_data.field_number,processed_data.cell_number]);

csvwrite(fullfile(output_dir,'average_cell_speed.csv'),processed_data.average_cell_speed);
csvwrite(fullfile(output_dir,'median_cell_speed.csv'),processed_data.average_cell_speed);

csvwrite(fullfile(output_dir,'degradation_area.csv'),processed_data.degradation_area);

csvwrite(fullfile(output_dir,'live_cells_per_image.csv'),processed_data.live_cells_per_image);

csvwrite(fullfile(output_dir,'corrected_final_diffs.csv'),processed_data.corrected_final_diff);

csvwrite(fullfile(output_dir,'gel_minus_surrounding.csv'),processed_data.gel_minus_surrounding);

degrade_indexes = find(processed_data.ever_degrade);
if (not(exist(fullfile(output_dir,'degrader'),'dir')))
    mkdir(fullfile(output_dir,'degrader'));
end

csvwrite(fullfile(output_dir,'degrader','degradation_area.csv'), ...
    processed_data.degradation_area(degrade_indexes));

csvwrite(fullfile(output_dir,'degrader','median_cell_speed.csv'), ...
    processed_data.median_cell_speed(degrade_indexes));

csvwrite(fullfile(output_dir,'degrader','corrected_final_diffs.csv'), ...
    processed_data.corrected_final_diff(degrade_indexes));

csvwrite(fullfile(output_dir,'degrader','fields_cells.csv'),...
     [processed_data.field_number(degrade_indexes),processed_data.cell_number(degrade_indexes)]);

function processed_data = process_raw_data(raw_data,varargin)

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

processed_data = struct();
if (isempty(strcmp('filter_set',i_p.UsingDefaults)))
    these_names = fieldnames(raw_data);
    for j=1:size(these_names,1)
        raw_data.(these_names{j}) = raw_data.(these_names{j})(i_p.Results.filter_set,:);
        processed_data.(these_names{j}) = raw_data.(these_names{j});
    end
end

processed_data.average_cell_speed = nanmean(processed_data.cell_speed,2);
processed_data.median_cell_speed = nanmedian(processed_data.cell_speed,2);

processed_data.live_cells = raw_data.tracking > 0;
processed_data.live_cells_per_image = sum(processed_data.live_cells);

processed_data.ever_degrade = [];
for i=1:size(raw_data.active_degrade,1)
    processed_data.ever_degrade = [processed_data.ever_degrade, any(raw_data.active_degrade(i,:))];
end

processed_data.has_degraded = zeros(size(raw_data.active_degrade));
for i=1:size(raw_data.active_degrade,1)
    for j = 1:size(raw_data.active_degrade,2)
        processed_data.has_degraded(i,j) = raw_data.active_degrade(i,j) | any(processed_data.has_degraded(i,1:j));
    end
end

processed_data.degraded_count = sum(processed_data.has_degraded);
processed_data.degrade_percentage = sum(processed_data.has_degraded)/size(raw_data.active_degrade,1);