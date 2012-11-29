function make_tracking_visualization(base_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('field_filter',0,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

addpath(genpath('../find_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

if (not(any(strcmp('field_filter',i_p.UsingDefaults))))
    fields = fields(i_p.Results.field_filter);
end

for field_num=1:length(fields)
    exp_dir = fullfile(base_dir,fields(field_num).name);
    image_dir = fullfile(exp_dir,'individual_pictures');
    single_image_dirs = dir(image_dir);

    single_image_dirs = single_image_dirs(3:end);
    
    tracking_file = fullfile(image_dir, single_image_dirs(1).name,filenames.tracking);
    if (not(exist(tracking_file,'file')))
        tracking_seq = zeros(1,size(single_image_dirs,1));
    else
        tracking_seq = csvread(tracking_file);
    end
    
    lineage_cmap = jet(size(tracking_seq,1));
    
    convert_avail = not(system('which convert > /dev/null'));
    for i_num = 1:length(single_image_dirs)
        current_data = read_in_file_set(fullfile(image_dir,single_image_dirs(i_num).name),filenames);
        
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
        
        output_file = fullfile(image_dir,single_image_dirs(i_num).name,filenames.tracking_vis);
        
        imwrite(highlighted_puncta,output_file);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Image Labeling
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if (not(convert_avail) || not(exist(tracking_file,'file'))), continue; end
        props = regionprops(current_data.labeled_cells,'Centroid');
        
        i_size = size(highlighted_puncta);
        for cell_num=1:length(props)
            if (props(cell_num).Centroid(1) > i_size(2)*0.95)
                props(cell_num).Centroid(1) = i_size(2)*0.95;
            end
            if (props(cell_num).Centroid(2) > i_size(1)*0.95)
                props(cell_num).Centroid(2) = i_size(1)*0.95;
            end
        end
        
        all_annotate = '';
        for cell_num = 1:length(ad_nums)
            pos_str = ['+',num2str(props(ad_nums(cell_num)).Centroid(1)),'+',num2str(props(ad_nums(cell_num)).Centroid(2))];
            label_str = [' "',num2str(tracking_nums(cell_num)),'"'];
            all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
        end
        command_str = ['convert ', output_file, ' -font VeraBd.ttf -pointsize 24 -fill ''rgb(255,0,0)''', ...
            all_annotate, ' ', output_file, '; '];
        system(command_str);
        1;
    end
    
    disp(['Done with ', exp_dir]);
end
toc;
