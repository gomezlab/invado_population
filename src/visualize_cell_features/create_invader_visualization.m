function create_invader_visualization(base_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('field_filter',0,@(x)isnumeric(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('..'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
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
    
    %toss out the '.' and '..' entries
    single_image_dirs = single_image_dirs(3:end);
    
    tracking_file = fullfile(image_dir, single_image_dirs(1).name,filenames.tracking);
    %check for the existance of a tracking file, if absent, there weren't any
    %cells in this field
    if (not(exist(tracking_file,'file')))
        disp('No tracking matrix found, assuming no cells in field');
        tracking_mat = zeros(size(single_image_dirs,1),1);
        active_degrade = zeros(size(single_image_dirs,1),1);
    else
        tracking_mat = csvread(tracking_file);
        active_degrade = csvread(fullfile(exp_dir,'cell_props','active_degrade.csv'));
        
        data_sets_to_read = {'Area','Gel_diff_minus_surrounding','Centroid_x','Centroid_y'};
        raw_data = struct();
        
        for i = 1:length(data_sets_to_read)
            data_dir = fullfile(image_dir, single_image_dirs(1).name,filenames.lineage_dir);
            raw_data.(data_sets_to_read{i}) = csvread(fullfile(data_dir,[data_sets_to_read{i}, '.csv']));
        end
        
        raw_data.corrected_final_gel_diffs = ...
            csvread(fullfile(image_dir,single_image_dirs(1).name,filenames.corrected_final_gel_diffs));
        
        longevity = sum(not(isnan(raw_data.Area)),2)/2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Building Visualization Images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i_num = 1:size(single_image_dirs,1)
        output_file = fullfile(image_dir,single_image_dirs(i_num).name, filenames.invader_vis);
        
        current_dir = fullfile(image_dir,single_image_dirs(i_num).name);
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
            short_lived = longevity(this_tracking_row) < 10;
            if (short_lived)
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
    convert_avail = not(system('which convert > /dev/null'));
    %if convert or the tracking matrix file aren't present, we either can't do
    %the labeling or we won't have any labels to add. Either way, time to
    %jump to the next image set.
    if (not(convert_avail) || not(exist(tracking_file,'file')))
        disp(['Done with ', exp_dir]);
        continue;
    end
    
    for i_num = 1:size(single_image_dirs,1)
        output_file = fullfile(image_dir,single_image_dirs(i_num).name, filenames.invader_vis);
        
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
        
        area = filtered_data.Area;
        gel_diff_percent = filtered_data.Gel_diff_minus_surrounding;
        
        all_annotate = '';
        for cell_num = 1:length(gel_diff_percent)
            cell_id = find(tracking_col == cell_num);
            pos_str = [' +',num2str(centroid(cell_num,1)),'+',num2str(centroid(cell_num,2))];
            top_line = sprintf('%d/%d', cell_id, area(cell_num));
            
            label_str = [' "', top_line,' \n', ...
                sprintf('%.2f%',gel_diff_percent(cell_num)),' \n', ...
                sprintf('%.2f%',raw_data.corrected_final_gel_diffs(cell_id)),'"'];
            all_annotate = [all_annotate, ' -annotate ', pos_str, label_str]; %#ok<AGROW>
        end
        command_str = ['convert ', output_file, ' -undercolor ''rgba(1,1,1,0.75)'' -font VeraBd.ttf -pointsize 16 -fill ''rgba(255,255,255,0.5)''', ...
            all_annotate, ' ', output_file, ';'];
        
        system(command_str);
    end
    
    disp(['Done with ', exp_dir]);
end

toc;