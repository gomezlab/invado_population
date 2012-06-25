function flat_field_correct_images(exp_dir, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath('matlab_scripts');
filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(exp_dir);
fields = filter_to_time_series(fields);

for j=1:length(fields)
    base_dir = fullfile(exp_dir,fields(j).name,'individual_pictures');
    
    image_dirs = dir(base_dir);
    
    assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
    assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
    assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>
    
    image_dirs = image_dirs(3:end);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Puncta Correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    puncta_flat_field = double(imread(fullfile(base_dir,image_dirs(1).name,filenames.puncta_median))); 
    for i = 1:size(image_dirs)
        if (exist(fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.puncta]),'file'))
            continue;
        end
        puncta_image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.puncta)));
        
        puncta_image = puncta_image - puncta_flat_field;
        puncta_image(puncta_image < 0) = 0;
        
%         movefile(fullfile(base_dir,image_dirs(i).name,filenames.puncta), ...
%             fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.puncta]));
        
        imwrite(uint16(puncta_image),fullfile(base_dir,image_dirs(i).name,filenames.puncta),'bitdepth',16);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gel Correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    gel_flat_field = double(imread(fullfile(base_dir,image_dirs(1).name,filenames.gel_median)));
    gel_flat_mean = mean(gel_flat_field(:));
    gel_flat_field = gel_flat_field - mean(gel_flat_field(:));
    gel_flat_field = gel_flat_field - (gel_flat_mean - 650);
    for i = 1:size(image_dirs)
        if (exist(fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.gel]),'file'))
            continue;
        end
        gel_image = double(imread(fullfile(base_dir,image_dirs(i).name,filenames.gel)));
        
        gel_image = gel_image - gel_flat_field;
        gel_image(gel_image < 0) = 0;
        
%         movefile(fullfile(base_dir,image_dirs(i).name,filenames.gel), ...
%             fullfile(base_dir,image_dirs(i).name,['uncorrected_', filenames.gel]));
        
        imwrite(uint16(gel_image),fullfile(base_dir,image_dirs(i).name,filenames.gel),'bitdepth',16);        
    end
    
    disp(['Done with ', fields(j).name]);
end
toc;
