function make_tracking_visualization(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

filenames = add_filenames_to_struct(struct());

tracking_seq = csvread(fullfile(base_dir,image_dirs(1).name,filenames.tracking)) + 1;
if (isempty(tracking_seq))
    tracking_seq = zeros(1,max_image_num);
end

lineage_cmap = jet(size(tracking_seq,1));

convert_avail = not(system('which convert'));
for i_num = 1:length(image_dirs)
    current_data = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);

    %thicken the cell perimeter, easier for visualization
    labeled_cells_perim_thick = thicken_perimeter(current_data.labeled_cells_perim,current_data.labeled_cells);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Image Creation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    live_rows = tracking_seq(:,i_num) > 0;
    ad_nums = tracking_seq(live_rows,i_num);
    tracking_nums = find(live_rows);
    
    %Build the unique lineage highlighted image
    this_cmap(ad_nums,:) = lineage_cmap(live_rows,:); %#ok<AGROW>    
    highlighted_puncta = create_highlighted_image(current_data.puncta_image_norm, ...
        labeled_cells_perim_thick,'color_map',this_cmap);
    
    output_file = fullfile(base_dir,image_dirs(i_num).name,filenames.tracking_vis);
    
    imwrite(highlighted_puncta,output_file);
    
	if (not(convert_avail)), continue; end
    props = regionprops(current_data.labeled_cells,'Centroid');
    
    i_size = size(highlighted_puncta);
    for i=1:length(props)
        if (props(i).Centroid(1) > i_size(2)*0.95)
            props(i).Centroid(1) = i_size(2)*0.95;
        end
        if (props(i).Centroid(2) > i_size(1)*0.95)
            props(i).Centroid(2) = i_size(1)*0.95;
        end
    end
    
    all_annotate = '';
    for i = 1:length(ad_nums)
        pos_str = ['+',num2str(props(ad_nums(i)).Centroid(1)),'+',num2str(props(ad_nums(i)).Centroid(2))];
        label_str = [' "',num2str(tracking_nums(i)),'"'];
        all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
    end
    command_str = ['convert ', output_file, ' -font VeraBd.ttf -pointsize 24 -fill ''rgb(255,0,0)''', ...
        all_annotate, ' ', output_file, '; '];
    system(command_str);
    1;
end
toc;
