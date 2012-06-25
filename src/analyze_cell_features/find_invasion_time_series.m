function find_invasion_time_series(results_dir,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i_p = inputParser;

i_p.addRequired('results_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(results_dir,varargin{:});

if (i_p.Results.debug == 1), profile off; profile on; end

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in the data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

exp_dirs = dir(fullfile(results_dir,'*'));
exp_dirs = filter_to_time_series(exp_dirs);

raw_data = struct();
for i=1:length(exp_dirs)
    files = struct();
    files.p_vals = fullfile(results_dir,exp_dirs(i).name,'adhesion_props','lin_time_series','Cell_gel_diff_p_val.csv');
    files.cell_diffs = fullfile(results_dir,exp_dirs(i).name,'adhesion_props','lin_time_series','Cell_gel_diff.csv');
    files.tracking = fullfile(results_dir,exp_dirs(i).name,'tracking_matrices','tracking_seq.csv');
    files.overlap_area = fullfile(results_dir,exp_dirs(i).name,'adhesion_props','lin_time_series','Overlap_region_size.csv');
    
    these_types = fieldnames(files);
    for j = 1:length(these_types)
        this_file = files.(these_types{j});
        
        %matlab doesn't like you to reference fields that haven't been
        %created, so create files that aren't present yet before loading
        %data in
        if(isempty(strmatch(these_types{j},fieldnames(raw_data))))
            raw_data.(these_types{j}) = [];
        end
        
        if (exist(this_file,'file'))
            temp = load(this_file);
            raw_data.(these_types{j}) = [raw_data.(these_types{j}); temp];
        else
            disp(['Can''t find ',this_file])
        end
    end
    
    %check that all the raw data files are the same size on each loop
    these_names = fieldnames(raw_data);
    poss_name_combinations = combnk(1:length(these_names),2);    
    for j=1:size(poss_name_combinations,1)
        assert(all(size(raw_data.(these_names{poss_name_combinations(j,1)})) == ...
                   size(raw_data.(these_names{poss_name_combinations(j,2)}))))
    end
end

%trim the first column off of all the raw data sets, this column is from
%the first image where degradation couldn't have been detected yet
these_names = fieldnames(raw_data);
for j=1:size(these_names,1)
    raw_data.(these_names{j}) = raw_data.(these_names{j})(:,2:end);
end

raw_data.cleaned_cell_diffs = raw_data.cell_diffs;
raw_data.cleaned_cell_diffs(isnan(raw_data.cleaned_cell_diffs)) = 0;
raw_data.cleaned_overlap_area = raw_data.overlap_area;
raw_data.cleaned_overlap_area(isnan(raw_data.cleaned_overlap_area)) = 0;

cell_diff_prods = raw_data.cleaned_cell_diffs.*raw_data.cleaned_overlap_area;

cleaned_sum = cumsum(cell_diff_prods,2);

cell_diff_prods_with_nan = cell_diff_prods;
cell_diff_prods_with_nan(isnan(raw_data.cleaned_cell_diffs)) = NaN;

correlations = double(corr(cell_diff_prods_with_nan'));
for i = 1:size(cell_diff_prods,1)
    correlations(i,i) = NaN;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process data files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
process_data = process_raw_data(raw_data);

%filtering based on longevity
min_hours = 10;
longev_filt = process_data.longevities >= min_hours;
longev_filtered_data = process_raw_data(raw_data,'filter_set',longev_filt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image Building
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_dir = fullfile(results_dir,'overall_results');

if (not(exist(output_dir,'dir')))
    mkdir(output_dir)
end

hist(process_data.longevities)
ylabel('Number of Objects')
xlabel('Lifetime of Objects')
title('NS Test 20X Tiling Longevities')
exportfig(gcf,fullfile(output_dir,'longevities.pdf'))
close gcf;

subplot(1,2,1);
time_points = linspace(0.5,length(process_data.degrade_percentage)*0.5,length(process_data.degrade_percentage));
plot(time_points,process_data.degrade_percentage);
hold on;
plot(time_points,longev_filtered_data.degrade_percentage,'r');
legend('No Filtering',['>=', num2str(min_hours), ' Hr Lifetime'],'Location','SouthEast')
xlabel('Time (hours)');
ylabel('Percentage of Cells that Have Degraded');

hold off;

subplot(1,2,2);
plot(log(time_points),process_data.degrade_percentage);
hold on;
plot(log(time_points),longev_filtered_data.degrade_percentage,'r');
xlabel('Time (Log Hours)');
ylabel('Percentage of Cells that Have Degraded');

exportfig(gcf,fullfile(output_dir,'test.pdf'),'Color','rgb','width',6,'height',3)
close gcf;

csvwrite(fullfile(output_dir,'degrade_percentage.csv'),longev_filtered_data.degrade_percentage);

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

if (isempty(strmatch('filter_set',i_p.UsingDefaults)))
    these_names = fieldnames(raw_data);
    for j=1:size(these_names,1)
        raw_data.(these_names{j}) = raw_data.(these_names{j})(i_p.Results.filter_set,:);
    end
end

bonferroni_correction = sum(sum(not(isnan(raw_data.p_vals))));

process_data = struct();

% process_data.active_degrade = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05/bonferroni_correction ...
%     & not(isnan(raw_data.cell_diffs)) & raw_data.cell_diffs < 0;

temp_cell_diffs = raw_data.cell_diffs(:);
temp_cell_diffs = temp_cell_diffs(not(isnan(temp_cell_diffs)));

process_data.active_degrade = not(isnan(raw_data.p_vals)) & raw_data.p_vals < 0.05/bonferroni_correction ...
    & not(isnan(raw_data.cell_diffs)) & raw_data.cell_diffs < 0;

process_data.live_cells = raw_data.tracking > -1;
process_data.longevities = sum(process_data.live_cells,2)/2;

process_data.ever_degrade = [];
for i=1:size(raw_data.tracking,1)
    process_data.ever_degrade = [process_data.ever_degrade, any(process_data.active_degrade(i,:))];
end

process_data.has_degraded = zeros(size(raw_data.tracking));
for i=1:size(raw_data.tracking,1)
    for j = 1:size(raw_data.tracking,2)
        process_data.has_degraded(i,j) = process_data.active_degrade(i,j) | any(process_data.has_degraded(i,1:j));
    end
end

process_data.degrade_percentage = sum(process_data.has_degraded)/size(raw_data.tracking,1);