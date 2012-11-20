function build_overall_data_set(exp_dir,varargin)
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
raw_data = struct('tracking',[],'field_number',[],'cell_number',[]);

prop_type_files = struct('active_degrade','active_degrade.csv', ...
    'longevity','longevity.csv', ...
    'degradation_area','degradation_area.csv',...
    'degradation_overall_rate','degradation_overall_rate.csv',...
    'degradation_rate_over_four','degradation_rate_over_four.csv',...
    'corrected_final_diff','corrected_final_gel_diffs.csv');
prop_names = fieldnames(prop_type_files);

time_series_files = struct('gel_diff_minus_surround','Gel_diff_minus_surrounding.csv',...
    'cell_speed','Cell_speed.csv');
ts_prop_names = fieldnames(time_series_files);

for j=1:length(fields)
    props_base = fullfile(exp_dir,fields(j).name,'cell_props');

    tracking_file = fullfile(exp_dir,fields(j).name,'tracking_matrices','tracking_seq.csv');
    if (not(exist(tracking_file,'file')))
        fprintf('No tracking file found in field %s, moving to next set.\n',fields(j).name);
        continue;
    end
    
    for prop=prop_names'
        prop = char(prop);
        if (isfield(raw_data,prop))
            raw_data.(prop) = [raw_data.(prop); csvread(fullfile(props_base,prop_type_files.(prop)))];
        else
            raw_data.(prop) = csvread(fullfile(props_base,prop_type_files.(prop)));
        end
    end
    
    for prop=ts_prop_names'
        prop = char(prop);
        if (isfield(raw_data,prop))
            raw_data.(prop) = [raw_data.(prop); csvread(fullfile(props_base,'lin_time_series',time_series_files.(prop)))];
        else
            raw_data.(prop) = csvread(fullfile(props_base,'lin_time_series',time_series_files.(prop)));
        end
    end
    
    tracking_data = csvread(tracking_file);
    raw_data.tracking = [raw_data.tracking; tracking_data];
    
    raw_data.field_number = [raw_data.field_number; j*ones(size(tracking_data,1),1)];
    raw_data.cell_number = [raw_data.cell_number; (1:size(tracking_data,1))'];
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All property file output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_file = fullfile(output_dir,'single_cell_properties.csv');
header_vals = {'field_number','cell_number','percent_matrix_degraded', ...
    'average_cell_speed','median_cell_speed','degradation_area', ...
    'degradation_overall_rate','degradation_rate_over_four','degrader'};

fid = fopen(out_file,'wt');
for i = 1:length(header_vals)
    if (i ~= length(header_vals))
        fprintf(fid,'%s,',header_vals{i});
    else
        fprintf(fid,'%s\n',header_vals{i});
    end
end
fclose(fid);

property_matrix = [processed_data.field_number,processed_data.cell_number, ...
    processed_data.corrected_final_diff,processed_data.average_cell_speed, ...
    processed_data.median_cell_speed,processed_data.degradation_area,...
    processed_data.degradation_overall_rate,processed_data.degradation_rate_over_four, ...
    processed_data.ever_degrade'];

dlmwrite(out_file,property_matrix,'-append');

csvwrite(fullfile(output_dir,'live_cells_per_image.csv'),processed_data.live_cells_per_image);
csvwrite(fullfile(output_dir,'gel_minus_surrounding.csv'),processed_data.gel_diff_minus_surround);
csvwrite(fullfile(output_dir,'corrected_final_diff.csv'),processed_data.corrected_final_diff);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output of Degrader Properties
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (not(exist(fullfile(output_dir,'degrader'),'dir')))
    mkdir(fullfile(output_dir,'degrader'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%All property file output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_file = fullfile(output_dir,'degrader','single_cell_properties.csv');

fid = fopen(out_file,'wt');
for i = 1:length(header_vals)
    if (i ~= length(header_vals))
        fprintf(fid,'%s,',header_vals{i});
    else
        fprintf(fid,'%s\n',header_vals{i});
    end
end
fclose(fid);

degrader_property_matrix = property_matrix(logical(processed_data.ever_degrade),:);

dlmwrite(out_file,degrader_property_matrix,'-append');


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
1;

processed_data.has_degraded = zeros(size(raw_data.active_degrade));
for i=1:size(raw_data.active_degrade,1)
    for j = 1:size(raw_data.active_degrade,2)
        processed_data.has_degraded(i,j) = raw_data.active_degrade(i,j) | any(processed_data.has_degraded(i,1:j));
    end
end

processed_data.degraded_count = sum(processed_data.has_degraded);
processed_data.degrade_percentage = sum(processed_data.has_degraded)/size(raw_data.active_degrade,1);