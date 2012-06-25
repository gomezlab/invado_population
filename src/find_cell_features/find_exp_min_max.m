function find_exp_min_max(base_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

gel_pix_vals = get_all_image_array(base_dir,filenames.gel);
gel_ranges = find_image_ranges(gel_pix_vals);

clear gel_pix_vals;

puncta_pix_vals = get_all_image_array(base_dir,filenames.puncta);
puncta_ranges = find_image_ranges(puncta_pix_vals);

clear gel_pix_vals;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Range Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);

gel_output_file = fullfile(output_base,image_nums(1).name, filenames.gel_range);
puncta_output_file = fullfile(output_base,image_nums(1).name, filenames.puncta_range);

[output_path] = fileparts(gel_output_file);

if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

csvwrite(gel_output_file,gel_ranges);
csvwrite(puncta_output_file,puncta_ranges);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ranges = find_image_ranges(pixel_vals)

pixel_vals = reshape(pixel_vals,[],1);
pixel_vals = sort(pixel_vals);

one_percent_data = round(size(pixel_vals,1)*0.01);

ranges = [pixel_vals(1),pixel_vals(end)];
for i=[0.5,1:5]
    ranges = cat(1,ranges, ...
        [pixel_vals(round(one_percent_data*i)), ...
         pixel_vals(end-round(one_percent_data*i))]);
end