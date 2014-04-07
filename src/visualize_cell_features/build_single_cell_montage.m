function build_single_cell_montage(exp_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('pixel_size',0.323,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);


i_p.parse(exp_dir,varargin{:});

addpath(genpath('../find_cell_features'));
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

try
    tracking_mat = csvread(fullfile(base_dir,image_dirs(1).name,filenames.tracking));
catch %#ok<CTCH>
    disp('Empty or missing tracking sequence, quiting.')
    return;
end
if (isempty(tracking_mat))
    tracking_mat = zeros(1,max_image_num);
end

active_degrade = csvread(fullfile(base_dir,image_dirs(1).name,filenames.active_degrade));

data_sets = cell(size(tracking_mat,2),1);
for i_num = 1:size(tracking_mat,2)
    data_sets{i_num} = read_in_file_set(fullfile(base_dir,image_dirs(i_num).name),filenames);
end

for cell_num = 1:size(tracking_mat,1)
    track_row = tracking_mat(cell_num,:);
    degrade_row = active_degrade(cell_num,:);
    
    %Don't make montages for cells that live less than 10 frames
    if (sum(track_row > 0) < 10)
        continue;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Determine extent of cell coverage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    cell_coverage = [];
    for i_num = 1:length(track_row)
        if (track_row(i_num) <= 0), continue; end
        
        this_cell = data_sets{i_num}.labeled_cells == track_row(i_num);
        if (size(cell_coverage,1) == 0)
            cell_coverage = this_cell;
        else
            cell_coverage = cell_coverage | this_cell;
        end
    end
    
    row_min_max = [find(sum(cell_coverage,2),1,'first') - 10, ...
        find(sum(cell_coverage,2),1,'last') + 10];
    col_min_max = [find(sum(cell_coverage),1,'first') - 10,...
        find(sum(cell_coverage),1,'last') + 10];
    
    if (row_min_max(1) < 1), row_min_max(1) = 1; end
    if (col_min_max(1) < 1), col_min_max(1) = 1; end
    
    if (row_min_max(2) > size(cell_coverage,1)), row_min_max(2) = size(cell_coverage,1); end
    if (col_min_max(2) > size(cell_coverage,2)), col_min_max(2) = size(cell_coverage,2); end
    
    data_sets_trimmed = cell(size(tracking_mat,2),1);
    for i_num = 1:length(track_row)
%         if (track_row(i_num) <= 0), continue; end

        data_sets_trimmed{i_num} = trim_all_images(data_sets{i_num},row_min_max,col_min_max);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Build visualization frames
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_frames = [];
    puncta_frames = [];
    
    %upper corners will hold the x position for image number labeling
    upper_corners = [];
    i_num_sequence = find(track_row);
    for i_num = 1:length(track_row)
        if (track_row(i_num) <= 0), continue; end

        if (size(upper_corners,1) == 0)
            upper_corners = 10;
        else
            upper_corners = [upper_corners, upper_corners(end) + 1 + size(gel_frames{1},2)]; %#ok<AGROW>
        end
                
        this_cell = data_sets_trimmed{i_num}.labeled_cells == track_row(i_num);
        this_cell_perim = data_sets_trimmed{i_num}.labeled_cells_perim == track_row(i_num);
        not_this_cell_perim = data_sets_trimmed{i_num}.labeled_cells_perim ~= track_row(i_num) & ...
            data_sets_trimmed{i_num}.labeled_cells_perim > 0;
        
        thick_perim = thicken_perimeter(this_cell_perim,this_cell);
        
        perim_color = [1,0,0];
        if (degrade_row(i_num) == 1)
            perim_color = [0,1,0];
        end
        
        gel_frame = create_highlighted_image(data_sets_trimmed{i_num}.gel_image_norm,thick_perim,'color_map',perim_color);
        gel_frame = create_highlighted_image(gel_frame,not_this_cell_perim,'color_map',[0,0,0.5]);
        gel_frames{length(gel_frames) + 1} = gel_frame; %#ok<AGROW>
        
        puncta_frame = create_highlighted_image(data_sets_trimmed{i_num}.puncta_image_norm,thick_perim,'color_map',perim_color);
        puncta_frame = create_highlighted_image(puncta_frame,not_this_cell_perim,'color_map',[0,0,0.5]);
        puncta_frames{length(puncta_frames) + 1} = puncta_frame; %#ok<AGROW>
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create final montage image and output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_montage = create_montage_image_set(gel_frames,'num_cols',length(gel_frames));
    puncta_montage = create_montage_image_set(puncta_frames,'num_cols',length(puncta_frames));
    spacer = 0.5*ones(1,size(gel_montage,2),3);
    
    full_montage = cat(1,puncta_montage,spacer,gel_montage);
    
    image_num = sprintf('%03d',cell_num);
    [~,field_dir] = fileparts(fileparts(base_dir));
    [~,exp_type] = fileparts(fileparts(fileparts(base_dir)));
    [~,date] = fileparts(fileparts(fileparts(fileparts(base_dir))));
    output_file = fullfile(base_dir,image_dirs(1).name,filenames.single_cell_dir,...
        [date,'-',exp_type(1:2),'-',field_dir,'-',image_num,'.png']);
%     output_file = fullfile(base_dir,image_dirs(1).name,filenames.single_cell_dir,[image_num,'.png']);
    
    if (not(exist(fileparts(output_file),'dir')))
        mkdir(fileparts(output_file))
    end
    
    imwrite(full_montage,output_file);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Add Image Number Labels
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    all_annotate = '';
    for i = 1:length(i_num_sequence)
        pos_str = ['+',num2str(upper_corners(i)),'+25'];
        label_str = [' "',num2str(i_num_sequence(i)),'"'];
        all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
    end
    command_str = ['convert ', output_file, ' -font VeraBd.ttf -pointsize 24 -fill ''rgba(255,255,255,0.5)''', ...
        all_annotate, ' ', output_file, '; '];
    system(command_str);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Sampled Montage
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_temp = gel_frames(round(linspace(1,length(gel_frames),10)));
    gel_montage = create_montage_image_set(gel_temp,'num_cols',length(gel_temp));
    
    puncta_temp = puncta_frames(round(linspace(1,length(puncta_frames),10)));
    puncta_montage = create_montage_image_set(puncta_temp,'num_cols',length(puncta_temp));
    spacer = 0.5*ones(1,size(gel_montage,2),3);
    
    full_montage = cat(1,puncta_montage,spacer,gel_montage);
    
    image_num = sprintf('%03d',cell_num);
    [~,field_dir] = fileparts(fileparts(base_dir));
    [~,exp_type] = fileparts(fileparts(fileparts(base_dir)));
    [~,date] = fileparts(fileparts(fileparts(fileparts(base_dir))));
    output_file = fullfile(base_dir,image_dirs(1).name,filenames.single_cell_dir,...
        [date,'-',exp_type(1:2),'-',field_dir,'-',image_num,'_sampled.png']);
%     output_file = fullfile(base_dir,image_dirs(1).name,filenames.single_cell_dir,[image_num,'.png']);
    
    if (not(exist(fileparts(output_file),'dir')))
        mkdir(fileparts(output_file))
    end
    
    imwrite(full_montage,output_file);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function images = trim_all_images(images,row_min_max,col_min_max)
image_size = size(images.gel_image);

names = fieldnames(images);
for i = 1:size(names)
    this_field = names{i};
    this_size = size(images.(this_field));
    %we only want to trim the images, not all the other misc properties
    %that come with the time point set of data, so check to make sure the
    %field is the size of the gel image and if so, trim
    if (length(this_size) == length(image_size) && ...
        all(this_size == image_size))
    
        images.(this_field) = images.(this_field)(row_min_max(1):row_min_max(2), ...
            col_min_max(1):col_min_max(2));
    end
end
