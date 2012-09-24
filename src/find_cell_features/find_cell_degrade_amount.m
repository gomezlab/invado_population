function find_cell_degrade_amount(base_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fields = dir(base_dir);
fields = filter_to_time_series(fields);

for i=1:length(fields)
    exp_dir = fullfile(base_dir,fields(i).name);
    image_dir = fullfile(exp_dir,'individual_pictures');
    
    image_dirs = dir(image_dir);
    
    assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
    assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
    assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
    
    image_dirs = image_dirs(3:end);
    
    tracking_file = fullfile(image_dir, image_dirs(1).name,filenames.tracking);
    %check for the existance of a tracking file, if absent, there weren't any
    %cells in this field, return from the function
    if (not(exist(tracking_file,'file')))
        disp('No tracking matrix found, assuming no cells in field');
        continue;
    else
        tracking_mat = csvread(tracking_file);
    end
    
    gel_range = csvread(fullfile(image_dir,image_dirs(1).name,filenames.gel_range));
    
    first_gel_image = double(imread(fullfile(image_dir,image_dirs(1).name,filenames.gel)));
    first_gel_image_trunc = first_gel_image;
    outside_range = first_gel_image > gel_range(2,2) | first_gel_image < gel_range(2,1);
    first_gel_image_trunc(outside_range) = NaN;
    
    
    final_gel_image = double(imread(fullfile(image_dir,image_dirs(end).name,filenames.gel)));
    final_gel_image_trunc = final_gel_image;
    outside_range = final_gel_image > gel_range(2,2) | final_gel_image < gel_range(2,1);
    final_gel_image_trunc(outside_range) = NaN;
    
    diff_image = final_gel_image_trunc - first_gel_image_trunc;
    
    no_cell_regions = imread(fullfile(image_dir,image_dirs(1).name,filenames.no_cells));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Count the number of times a cell was seen in each pixel
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cell_hit_counts = cell(size(tracking_mat,1),1);
    
    for i_num = 1:size(tracking_mat,2)
        current_dir = fullfile(image_dir,image_dirs(i_num).name);
        labeled_cells = imread(fullfile(current_dir,filenames.labeled_cell_mask));
        
        for cell_num = 1:size(tracking_mat,1)
            if (tracking_mat(cell_num,i_num) == 0), continue; end
            
            if (isempty(cell_hit_counts{cell_num}))
                cell_hit_counts{cell_num} = zeros(size(labeled_cells,1),size(labeled_cells,2));
            end
            
            cell_hit_counts{cell_num} = cell_hit_counts{cell_num} + double(labeled_cells == tracking_mat(cell_num,i_num));
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine what the matrix does underneath each cell
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    diffs = zeros(length(cell_hit_counts),1);
    corrected_diffs = zeros(length(cell_hit_counts),1);
    
    for cell_num = 1:length(cell_hit_counts)
        cell_extent = cell_hit_counts{cell_num} > 10;
        
        surrounding_cell_extent = imdilate(cell_extent,strel('disk',40)) & not(cell_extent);
        
        diff_vals = diff_image(cell_extent);
        surrounding_diff_pixels = diff_image(surrounding_cell_extent);
        
        nan_count = sum(isnan(diff_vals));
        nan_surrounding_count = sum(isnan(surrounding_diff_pixels));
        
        temp = zeros(size(cell_extent));
        temp(cell_extent) = 1;
        temp(no_cell_regions) = 2;
        temp(surrounding_cell_extent) = 3;
        
        if (nan_count > length(diff_vals)*0.5 || ...
            nan_surrounding_count > length(surrounding_diff_pixels)*0.5)
            diffs(cell_num) = NaN;
            corrected_diffs(cell_num) = NaN;
        else
            first_intensity = nanmean(first_gel_image_trunc(cell_extent));
            diffs(cell_num) = 100*(nanmean(diff_vals)/first_intensity);
            corrected_diffs(cell_num) = 100*((nanmean(diff_vals) - nanmean(surrounding_diff_pixels))/first_intensity);
        end
    end
    
    csvwrite(fullfile(image_dir,image_dirs(i).name,filenames.final_gel_diffs),diffs);
    csvwrite(fullfile(image_dir,image_dirs(i).name,filenames.corrected_final_gel_diffs),corrected_diffs);
    disp(['Done with ', exp_dir]);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%