function build_raw_montage_file_set(base_dir,target_file,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('base_dir',@(x)exist(x,'dir') == 7);
i_p.addRequired('target_file',@(x)ischar(x)||isarray(x));
i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(base_dir,target_file,varargin{:});

addpath(genpath('..'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fields = dir(base_dir);
fields = filter_to_time_series(fields);

images_composite = cell(0);
base_image_size = [NaN, NaN];

total_images = length(dir(fullfile(base_dir,fields(1).name,'individual_pictures'))) - 2;
all_pix_vals = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

any_three_color = 0;
%build up composite images by scanning through the image number, so time
%step, then each field
for time = 1:total_images
    images = cell(0);
    
    for j=1:length(fields)
        image_base = fullfile(base_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(time).name, target_file);
        
        %deal with the case where the image does not exist
        if(not(exist(image_file,'file')))
            continue;
        end
        [images{j},map] = imread(image_file);
        %deal with indexed or single layer png images, occasionaly produced
        %by imagemagick
        if (size(map,1) > 0)
            images{j} = uint8(256*ind2rgb(images{j},map));
        end
        
        base_image_size = size(images{j});
        if (length(base_image_size) > 2)
            any_three_color = 1;
        end
    end
    
    %scan back through the images collected filling in missing files with
    %Nan blanks
    for j=1:length(fields)
        if(length(images) < j || isempty(images{j}))
            images{j} = NaN*ones(base_image_size);
        end
        if (any_three_color && length(size(images{j})) < 3)
            images{j} = cat(3,images{j},images{j},images{j});
        end
    end
    
    %build up composite image, cat first 5 images in one column
    images_composite{time} = cat(1,images{1:5});
    for col_num=2:(length(images)/5)
        image_nums = ((col_num-1)*5+1):(col_num*5);
        %if this is an even row, reverse the image order to compensate for
        %the zig-zag pattern of image collection
        if (mod(col_num,2) == 0)
            image_nums = image_nums(end:-1:1);
        end
        %then cat each subsequent five images in the next row
        images_composite{time} = cat(2,images_composite{time}, cat(1,images{image_nums}));
    end
    
    %finally resize to a reasonable image size
    images_composite{time} = imresize(images_composite{time},[800,NaN]);
    all_pix_vals = [all_pix_vals; images_composite{time}(:)]; %#ok<AGROW>
    
    if (mod(time,10) == 0)
        disp(['Done processing: ',num2str(time), ' of ', target_file])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalize Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

all_pix_vals = sort(all_pix_vals);
%we only need to do image normalization if the images are not already in
%color (i.e. have more than 2 dimensions)
if (ndims(images_composite{1}) < 3)
    %toss out the top 0.05% of pixels for the image min/max creation
    end_trim_amount = round(length(all_pix_vals)*0.0005);
    images_min_max = double([all_pix_vals(1), all_pix_vals(end-end_trim_amount)]);
    for i=1:length(images_composite)
        %scale and then truncate the composite images
        images_composite{i} = (double(images_composite{i}) - images_min_max(1))/range(images_min_max);
        images_composite{i}(images_composite{i} > images_min_max(2)) = images_min_max(2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output Composite Images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[pathstr, name] = fileparts(target_file); %#ok<ASGLU>
output_dir = fullfile(base_dir,'montage',name);

if (not(exist(fullfile(output_dir),'dir')))
    mkdir(fullfile(output_dir));
end

for i=1:length(images_composite)
    padding_length = num2str(length(num2str(total_images)));
    padded_i_num = sprintf(['%0',padding_length,'d'],i);
    
    imwrite(images_composite{i}, fullfile(output_dir,[padded_i_num,'.png']));
    if (mod(i,10) == 0)
        disp(['Done writing: ',num2str(i),' of ', target_file])
    end
end
toc;
