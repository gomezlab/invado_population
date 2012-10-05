function find_cell_mask_properties(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

data_points = 0;
all_tracking_props = cell(size(image_dirs,1),1);
all_cell_props = cell(size(image_dirs,1),1);

final_data = read_in_file_set(fullfile(base_dir,image_dirs(size(image_dirs,1)).name),filenames);

for i_num = 1:size(image_dirs,1)
    if (i_num > 1)
        prior_dir = fullfile(base_dir,image_dirs(i_num-1).name);
        if (exist('current_data','var'))
            prior_data = current_data;
        else
            prior_data = read_in_file_set(prior_dir,filenames);
        end
    else
        prior_dir = fullfile(base_dir,image_dirs(i_num).name);
        prior_data = read_in_file_set(prior_dir,filenames);
    end
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    gel_junk_threshold = csvread(fullfile(current_dir, filenames.gel_junk_threshold));        
    
    current_data.gel_junk = current_data.gel_image > gel_junk_threshold;
    prior_data.gel_junk = prior_data.gel_image > gel_junk_threshold;
    
    cell_props = collect_cell_properties(current_data,prior_data,final_data,'debug',i_p.Results.debug);
    all_cell_props{i_num} = cell_props;
    
    data_points = data_points + length(cell_props);
    
    if (i_num ~= 1)
        tracking_props = collect_tracking_properties(current_data,prior_data,'debug',i_p.Results.debug);
        
        all_tracking_props{i_num-1} = tracking_props;
    end
    
    if (i_num == size(image_dirs,1))
        tracking_props = collect_tracking_properties(current_data,current_data,'debug',i_p.Results.debug);
        all_tracking_props{i_num} = tracking_props;
    end
    
    %make a diagnostic figure showing where the cell masks overlap with the
    %next image
    temp = double(current_data.cell_mask);
    temp(prior_data.cell_mask) = 2;
    temp(prior_data.cell_mask & current_data.cell_mask) = 3;
    temp(current_data.gel_junk | prior_data.gel_junk) = 4;
    imwrite(label2rgb(temp),fullfile(current_dir,'cell_overlaps.png'));
    
    diagnostic_boundaries = get_obj_perims(temp);
    %     diagnostic_boundaries(current_data.gel_junk | prior_data.gel_junk) = 4;
    
    diag_cmap = [[0,1,0];[0,0,1];[0,1,1];[1,0,0]];
    %make an image showing where the cells are located overlayed on the gel
    %image
    gel_range_norm = create_highlighted_image(current_data.gel_image_norm,diagnostic_boundaries,'color_map',diag_cmap);
    imwrite(gel_range_norm,fullfile(current_dir,'gel_highlights.png'));
    
    if (mod(i_num,10)==0)
        disp(['Done processing image number: ',num2str(i_num)])
    end
end

save(fullfile(base_dir,image_dirs(1).name,filenames.tracking_raw),'all_tracking_props');
save(fullfile(base_dir,image_dirs(1).name,filenames.cell_props),'all_cell_props');
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function obj_perims = get_obj_perims(label_mat)

obj_perims = zeros(size(label_mat));

if (max(label_mat(:)) == 0)
    return;
end

for i = 1:max(label_mat(:))
    this_obj = label_mat == i;
    this_perim = bwperim(this_obj);
    obj_perims(this_perim) = i;
end

function cell_props = collect_cell_properties(current_data,prior_data,final_data,varargin)
% COLLECT_cell_PROPERTIES    using the identified cells, various
%                                properties are collected concerning the
%                                morphology and physical properties of the
%                                cells

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('current_data',@isstruct);
i_p.addRequired('prior_data',@isstruct);
i_p.addRequired('final_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prior_data,final_data,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cell_props = regionprops(current_data.labeled_cells,current_data.gel_image, ...
    'Area','Centroid','MeanIntensity','PixelValues');

for i=1:max(current_data.labeled_cells(:))
    cell_props(i).StdIntensity = std([cell_props(i).PixelValues]);
end
cell_props = rmfield(cell_props,'PixelValues');

if (isempty(cell_props))
    cell_props(1).Overlap_area = [];
    cell_props(1).Overlap_percent = [];
    
    cell_props(1).Cell_gel_before = [];
    cell_props(1).Cell_gel_after = [];
    
    cell_props(1).Cell_gel_diff = [];
    cell_props(1).Cell_gel_diff_p_val = [];
    cell_props(1).Cell_gel_diff_median = [];
    cell_props(1).Cell_gel_diff_total = [];
    cell_props(1).Cell_gel_diff_percent = [];
    
    cell_props(1).Surrounding_diff_percent = [];
    cell_props(1).Gel_diff_minus_surrounding = [];
else
    [cell_props.Overlap_area] = deal(NaN);
    [cell_props.Overlap_percent] = deal(NaN);
    
    [cell_props.Cell_gel_before] = deal(NaN);
    [cell_props.Cell_gel_after] = deal(NaN);
    
    [cell_props.Cell_gel_diff] = deal(NaN);
    [cell_props.Cell_gel_diff_p_val] = deal(NaN);
    [cell_props.Cell_gel_diff_median] = deal(NaN);
    [cell_props.Cell_gel_diff_total] = deal(NaN);
    [cell_props.Cell_gel_diff_percent] = deal(NaN);

    [cell_props.Surrounding_diff_percent] = deal(NaN);
    [cell_props.Gel_diff_minus_surrounding] = deal(NaN);
end

%when the first image is both the prior and current data, we only want the
%area and centroid variables to contain any data, the other properties are
%not relavent for the first image
if (all(current_data.gel_image(:) == prior_data.gel_image(:)) == 1)
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Properites Always Extracted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diff_image = current_data.gel_image - prior_data.gel_image;

for i=1:max(current_data.labeled_cells(:))
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_cell = current_data.labeled_cells;
    this_cell(current_data.labeled_cells ~= i) = 0;
    this_cell = logical(this_cell);
    
    prev_cells = prior_data.cell_mask;

    this_surrounding_region = imdilate(this_cell,strel('disk',20)) & not(this_cell) & not(prev_cells);    
    
    overlap_region = this_cell & prev_cells & not(current_data.gel_junk) & ...
        not(prior_data.gel_junk);
    cell_props(i).Overlap_area = sum(overlap_region(:));
    
    cell_props(i).Overlap_percent = cell_props(i).Overlap_area/cell_props(i).Area;
    %If the percentage of valid overlap pixels falls below 40%, we want to
    %leave all the other properties as defaults, so...
    if (cell_props(i).Overlap_percent < 0.4)
        continue;
    end
    
    %Cell overlap differences
    differences = diff_image(overlap_region);
    
    cell_props(i).Cell_gel_before = mean(prior_data.gel_image(overlap_region));
    cell_props(i).Cell_gel_after = mean(current_data.gel_image(overlap_region));
    
    cell_props(i).Cell_gel_diff = mean(differences);
    [~,p] = ttest(differences);
    cell_props(i).Cell_gel_diff_p_val = p;
    cell_props(i).Cell_gel_diff_median = median(differences);
    cell_props(i).Cell_gel_diff_total = sum(differences);
    
    cell_props(i).Cell_gel_diff_percent = 100*(cell_props(i).Cell_gel_diff/cell_props(i).Cell_gel_before);
    
    %Surrounding differences
    surrounding_diff = diff_image(this_surrounding_region);
    
    cell_props(i).Surrounding_diff_percent = 100*(mean(surrounding_diff)/cell_props(i).Cell_gel_before);
    
    cell_props(i).Gel_diff_minus_surrounding = cell_props(i).Cell_gel_diff_percent - ...
        cell_props(i).Surrounding_diff_percent;
    
    %single cell diagnostics
    if (i_p.Results.debug)
        c_extent = [find(sum(this_cell,2), 1 ),find(sum(this_cell,2), 1, 'last' ), ...
            find(sum(this_cell), 1 ),find(sum(this_cell), 1, 'last' )];
        c_extent(1:2:3) = c_extent(1:2:3) - 10;
        c_extent(2:2:4) = c_extent(2:2:4) + 10;
        c_extent(c_extent <= 0) = 1;
        if (c_extent(2) > size(this_cell,1)), c_extent(2) = size(this_cell,1); end
        if (c_extent(4) > size(this_cell,2)), c_extent(4) = size(this_cell,2); end
        c_extent = struct('row',(c_extent(1):c_extent(2)),'col',(c_extent(3):c_extent(4)));
        
        prior_puncta = normalize_image(prior_data.puncta_image(c_extent.row,c_extent.col));
        prior_gel = prior_data.gel_image(c_extent.row,c_extent.col)*prior_data.intensity_correction;
        
        current_puncta = normalize_image(current_data.puncta_image(c_extent.row,c_extent.col));
        current_puncta = (current_puncta - min(current_puncta(:)))/range(current_puncta(:));
        current_gel = current_data.gel_image(c_extent.row,c_extent.col)*current_data.intensity_correction;
        
        gel_min = min([prior_gel(:);current_gel(:)]);
        gel_max = max([prior_gel(:);current_gel(:)]);
        prior_gel = normalize_image(prior_gel,[gel_min,gel_max]);
        current_gel = normalize_image(current_gel,[gel_min,gel_max]);
        
        prior_both = cat(1,prior_puncta,prior_gel);
        current_both = cat(1,current_puncta,current_gel);
        
        diff_gel = normalize_image(current_gel - prior_gel);
        
        temp = double(this_cell);
        temp(prev_cells) = 2;
        temp(overlap_region) = 3;
        subplot(2,2,1); imshow(label2rgb(temp));
        subplot(2,2,2); hist(differences);
        subplot(2,2,3); imshow(diff_gel);
        1;
    end
end

function prior_props = collect_tracking_properties(current_data,prior_data,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('current_data',@isstruct);
i_p.addRequired('prior_data',@isstruct);

i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(current_data,prior_data,varargin{:});

current_props = regionprops(current_data.labeled_cells,'Centroid');
prior_props = regionprops(prior_data.labeled_cells,'Centroid','Area');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (all(current_data.gel_image(:) == prior_data.gel_image(:)) == 1)
    return;
end

for i=1:max(prior_data.labeled_cells(:))
    for j=1:max(current_data.labeled_cells(:))
        prior_props(i).Cent_dist(j) = sqrt((current_props(j).Centroid(1) - prior_props(i).Centroid(1))^2 + ...
            (current_props(j).Centroid(2) - prior_props(i).Centroid(2))^2);
    end
end

for i=1:max(prior_data.labeled_cells(:))
    %this bit of code isolates a single object as a logical image and
    %then builds another logical image of the region around the object,
    %excluding certain areas
    this_cell = prior_data.labeled_cells;
    this_cell(prior_data.labeled_cells ~= i) = 0;
    this_cell = logical(this_cell);
    
    this_cell_area = sum(sum(this_cell));
    next_overlaping_cells = current_data.labeled_cells(this_cell);
    if (max(current_data.labeled_cells(:)) == 0)
        continue;
    end
    
    prior_props(i).Pix_sim = zeros(1,max(current_data.labeled_cells(:)));
    overlap_nums = unique(next_overlaping_cells);
    if (overlap_nums(1) == 0)
        overlap_nums = overlap_nums(2:end);
    end
    for j=1:length(overlap_nums)
        this_overlap_num = overlap_nums(j);
        
        overlap_area = sum(next_overlaping_cells == this_overlap_num);
        prior_props(i).Pix_sim(this_overlap_num) = overlap_area/this_cell_area;
    end
end
