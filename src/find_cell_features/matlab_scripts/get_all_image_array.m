function images = get_all_image_array(exp_dir,file_name)
% This function gathers all the images of a specific type into a single
% array.

fields = dir(exp_dir);
fields = filter_to_time_series(fields);

total_images = length(dir(fullfile(exp_dir,fields(1).name,'individual_pictures'))) - 2;

images_cell = cell(0);
for i = 1:total_images
    for j=1:length(fields)
        image_base = fullfile(exp_dir,fields(j).name,'individual_pictures');
        image_nums = dir(image_base);
        image_nums = image_nums(3:end);
        image_file = fullfile(image_base,image_nums(i).name, file_name);
                
        images_cell{length(images_cell)+1} = imread(image_file);
    end
        
    if (mod(i,10) == 0)
        disp(['Done reading: ',num2str(i), '/', num2str(total_images), ' for ', file_name])
    end
end

images = zeros([size(images_cell{1}),length(images_cell)]);

for i = 1:length(images_cell)
    images(:,:,i) = images_cell{i};
end