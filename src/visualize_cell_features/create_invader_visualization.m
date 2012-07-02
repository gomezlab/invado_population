function create_invader_visualization(field_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('field_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(field_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
base_dir = fullfile(field_dir,'individual_pictures');
image_dirs = dir(base_dir);
assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
image_dirs = image_dirs(3:end);

tracking_file = fullfile(base_dir, image_dirs(1).name,filenames.tracking);
%check for the existance of a tracking file, if absent, there weren't any
%cells in this field
if (not(exist(tracking_file,'file')))
    disp('No tracking matrix found, assuming no cells in field');
    tracking_mat = zeros(size(image_dirs,1),1);
    active_degrade = zeros(size(image_dirs,1),1);
else
    tracking_mat = csvread(tracking_file);
    active_degrade = csvread(fullfile(field_dir,'cell_props','active_degrade.csv'));
    
    data_sets_to_read = {'Cell_gel_diff','Cell_gel_diff_median','Cell_gel_diff_p_val','Area' ...
        'Cell_gel_diff_total','Cell_gel_diff_percent','Centroid_x','Centroid_y','Overlap_area'};
    raw_data = struct();
    
    for i = 1:length(data_sets_to_read)
        data_dir = fullfile(base_dir, image_dirs(1).name,filenames.lineage_dir);
        raw_data.(data_sets_to_read{i}) = csvread(fullfile(data_dir,[data_sets_to_read{i}, '.csv']));
    end
    
    longevity = sum(not(isnan(raw_data.Area)),2)/2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Building Visualization Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i_num = 1:size(image_dirs,1)
    output_file = fullfile(base_dir,image_dirs(i_num).name, filenames.invader_vis);
    
    current_dir = fullfile(base_dir,image_dirs(i_num).name);
    current_data = read_in_file_set(current_dir,filenames);
    
    degrade_marked = zeros(size(current_data.labeled_cells));
    for cell_num=1:max(current_data.labeled_cells_perim(:))
        this_tracking_row = tracking_mat(:,i_num) == cell_num;
        assert(sum(this_tracking_row) == 1);
        
        degrade_status = active_degrade(this_tracking_row,i_num);

        
        %use 1 for non-degraders, 2 for degraders
        if (degrade_status == 0)
            degrade_marked(current_data.labeled_cells_perim == cell_num) = 1;
        elseif (degrade_status == 1)
            degrade_marked(current_data.labeled_cells_perim == cell_num) = 2;
        else
            disp('Found unexpected degrade status');
        end
        
        %use 3 for short lived cells, without regard for the degradation
        %status
        long_lived = longevity(this_tracking_row) >= 10;        
        if (not(long_lived))
            degrade_marked(current_data.labeled_cells_perim == cell_num) = 3;
        end
        
    end
    
    c_map = [[1,0,0];[0,1,0];[95/255,0,128/255]];
    
    thick_degrade_marked = thicken_perimeter(degrade_marked,current_data.labeled_cells);
    
    degrade_highlights = create_highlighted_image(current_data.gel_image_norm,thick_degrade_marked,'color_map',c_map);
    imwrite(degrade_highlights,output_file);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adding Degradation Annotations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
convert_avail = not(system('which convert'));
%if convert or the tracking matrix file aren't present, we either can't do
%the labeling or we won't have any labels to add. Either way, time to
%return.
if (not(convert_avail) || not(exist(tracking_file,'file')))
    return;
end

for i_num = 1:size(image_dirs,1)
    output_file = fullfile(base_dir,image_dirs(i_num).name, filenames.invader_vis);
    
    if (not(exist(output_file,'file')))
        continue;
    end
    
    img_size = size(imread(output_file));
    
    tracking_col = tracking_mat(:,i_num);
    
    filtered_data = struct();
    for i = 1:length(data_sets_to_read)
        temp = NaN(sum(tracking_col > 0),1);
        for j=1:length(tracking_col)
            if (tracking_col(j) < 1)
                continue;
            end
            temp(tracking_col(j)) = raw_data.(data_sets_to_read{i})(j, i_num);
        end
        filtered_data.(data_sets_to_read{i}) = temp;
    end
    
    filtered_data.Centroid_x(filtered_data.Centroid_x > 0.9*img_size(2)) = 0.9*img_size(2);
    filtered_data.Centroid_y(filtered_data.Centroid_y > 0.9*img_size(1)) = 0.9*img_size(1);
    filtered_data.Centroid_y(filtered_data.Centroid_y < 0.1*img_size(1)) = 0.1*img_size(1);
    centroid = [filtered_data.Centroid_x,filtered_data.Centroid_y];
    
    gel_diff = filtered_data.Cell_gel_diff;
    area = filtered_data.Area;
    gel_diff_percent = filtered_data.Cell_gel_diff_percent;
    gel_diff_median = filtered_data.Cell_gel_diff_median;
    
    all_annotate = '';
    for cell_num = 1:length(gel_diff)
        pos_str = [' +',num2str(centroid(cell_num,1)),'+',num2str(centroid(cell_num,2))];
        top_line = sprintf('%d/%d/%d',find(tracking_col == cell_num), ...
            area(cell_num),filtered_data.Overlap_area(cell_num));
        label_str = [' "', top_line,' \n', ...
            sprintf('%.2f',gel_diff_percent(cell_num)), '% \n', ...
            sprintf('%.2f',gel_diff(cell_num)),'"'];
        all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
    end
    command_str = ['convert ', output_file, ' -undercolor ''rgba(1,1,1,0.75)'' -font VeraBd.ttf -pointsize 16 -fill ''rgba(255,255,255,0.5)''', ...
        all_annotate, ' ', output_file, ';'];
    
    system(command_str);
end

toc;
