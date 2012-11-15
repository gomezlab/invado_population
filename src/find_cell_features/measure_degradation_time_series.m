function measure_degradation_time_series(base_dir,varargin)

tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('field_filter',0,@(x)isnumeric(x));
i_p.addParamValue('min_fraction_decrease',0.2,@(x)isnumeric(x));
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

if (not(any(strcmp('field_filter',i_p.UsingDefaults))))
   fields = fields(i_p.Results.field_filter); 
end

for i=1:length(fields)
    exp_dir = fullfile(base_dir,fields(i).name);
    image_dir = fullfile(exp_dir,'individual_pictures');
    
    single_image_folders = dir(image_dir);
    
    assert(strcmp(single_image_folders(1).name, '.'), 'Error: expected "." to be first string in the dir command')
    assert(strcmp(single_image_folders(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
    assert(str2num(single_image_folders(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
    
    single_image_folders = single_image_folders(3:end);
    
    tracking_file = fullfile(image_dir, single_image_folders(1).name,filenames.tracking);
    %check for the existance of a tracking file, if absent, there weren't any
    %cells in this field, return from the function
    if (not(exist(tracking_file,'file')))
        disp('No tracking matrix found, assuming no cells in field');
        continue;
    else
        tracking_mat = csvread(tracking_file);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reading and Processing first image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gel_junk_threshold = csvread(fullfile(image_dir,single_image_folders(1).name, filenames.gel_junk_threshold));        
    
    first_gel_image = double(imread(fullfile(image_dir,single_image_folders(1).name,filenames.gel)));
    first_gel_image_trunc = first_gel_image;
    junk_area = first_gel_image > gel_junk_threshold | ...
        first_gel_image < 30;
    first_gel_image_trunc(junk_area) = NaN;
    
    gel_junk_regions = logical(imread(fullfile(image_dir,single_image_folders(1).name,filenames.gel_junk)));
    
    degraded_areas = zeros(size(first_gel_image));
    degraded_votes = zeros(size(first_gel_image));
    degraded_time = zeros(size(first_gel_image));
    degraded_areas_labeled = zeros(size(first_gel_image));
    
    single_cell_degraded_area = zeros(size(tracking_mat));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Processing Remaining Images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    for i_num = 2:size(tracking_mat,2)
        current_dir = fullfile(image_dir,single_image_folders(i_num).name);
        data_set = read_in_file_set(current_dir,filenames);
        
        gel_image_trunc = data_set.gel_image;
        gel_image_trunc(gel_image_trunc > gel_junk_threshold) = NaN;
        
        %Processing the 
        degraded_area_so_far = gel_image_trunc < first_gel_image_trunc*(1-i_p.Results.min_fraction_decrease);
        degraded_area_so_far(gel_junk_regions) = 0;
        degraded_area_so_far(isnan(first_gel_image_trunc)) = 0;
        degraded_area_so_far(isnan(degraded_area_so_far)) = 0;
        degraded_votes = degraded_votes + degraded_area_so_far;
        
        new_degrade = not(degraded_areas) & degraded_area_so_far;
        degraded_time(new_degrade) = i_num;
        
        degraded_areas = degraded_areas | degraded_area_so_far;
        
        %Debugging images
%         gel_image_trunc_norm = data_set.gel_image_norm;
%         gel_image_trunc_norm(isnan(gel_image_trunc)) = 0;
%         temp = degraded_areas;
%         temp(new_degrade) = 2;
%         close;
%         imshow(create_highlighted_image(gel_image_trunc_norm,temp,'mix_percent',1));
%         pause(0.1);        

        for cell_num = 1:size(tracking_mat,1)
            %zeros in the tracking mat indicate no cell
            if (tracking_mat(cell_num,i_num) == 0)
                continue;
            end
            this_cell_degrade = new_degrade & data_set.labeled_cells == tracking_mat(cell_num,i_num);
            single_cell_degraded_area(cell_num,i_num) = sum(this_cell_degrade(:));
            degraded_areas_labeled(this_cell_degrade) = cell_num;
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Final Results Processing and Output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    cumulative_degraded_area = cumsum(single_cell_degraded_area,2);
    final_degrade_area = cumulative_degraded_area(:,end);
    
    % plot(cumsum(single_cell_degraded_area,2)');
    
    output_dir = fullfile(image_dir,single_image_folders(i).name,filenames.lineage_dir);
    csvwrite(fullfile(output_dir,'area_degraded.csv'),single_cell_degraded_area);
    csvwrite(fullfile(output_dir,'cumul_area_degraded.csv'),cumulative_degraded_area);
    
    csvwrite(fullfile(image_dir,single_image_folders(i).name,filenames.degrade_areas),final_degrade_area);
    
    disp(['Done with ', exp_dir]);
end
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%