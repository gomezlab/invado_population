function find_median_images(base_dir,varargin)
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

output_base = fullfile(base_dir,fields(1).name,'individual_pictures');
image_nums = dir(output_base);
image_nums = image_nums(3:end);

output_file = fullfile(output_base,image_nums(1).name, filenames.puncta_median);
[output_path] = fileparts(output_file);
if (not(exist(fullfile(output_path),'dir')))
    mkdir(fullfile(output_path));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Puncta Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

images_puncta = get_all_image_array(base_dir,filenames.puncta);
puncta_median = median(images_puncta,3);

imwrite(uint16(puncta_median),output_file,'bitdepth',16);

clear images_puncta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gel Processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

images_gel = get_all_image_array(base_dir,filenames.gel);
gel_median = median(images_gel,3);

output_file = fullfile(output_base,image_nums(1).name, filenames.gel_median);
imwrite(uint16(gel_median),output_file,'bitdepth',16);

clear images_gel;
toc;