function track_cells(exp_dir,varargin)
%CREATE_CELL_MASK_IMAGE   Gather and write the cell mask from a
%                         fluorescence image
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

%Add the folder with all the scripts used in this master program
addpath(genpath('matlab_scripts'));
% addpath(genpath('../visualize_cell_features'));

filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

base_dir = fullfile(exp_dir,'individual_pictures');

image_dirs = dir(base_dir);

assert(strcmp(image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

image_dirs = image_dirs(3:end);

load(fullfile(base_dir,image_dirs(1).name,filenames.tracking_raw));

%fill in the possible linked objects based on the number of values in the
%area
for i_num=1:length(all_tracking_props) %#ok<NODEF>
    area = [all_tracking_props{i_num}.Area];
    
    for i=1:size(area,2)
        all_tracking_props{i_num}(i).assigned = 0;
    end
    try
        pix_sim = reshape([all_tracking_props{i_num}.Pix_sim],[],size(area,2));
        for i=1:size(area,2)
            all_tracking_props{i_num}(i).next_obj = [];
        end
    catch
        continue;
    end
end

%filter the tracking props on the pixel similary measure
for i_num=1:(length(all_tracking_props)-1)
    pix_sim = reshape([all_tracking_props{i_num}.Pix_sim],[],length(all_tracking_props{i_num}))';
    
    high_pix_sim = pix_sim > 0.5;
    
    [start_cell_hits, end_cell_hits] = find(high_pix_sim);
    
    for i = 1:length(start_cell_hits)
        start_cell = start_cell_hits(i);
        end_cell = end_cell_hits(i);
        if (sum(start_cell == start_cell_hits) == 1 && ...
                sum(end_cell == end_cell_hits) == 1)            
            
            all_tracking_props{i_num}(start_cell).next_obj = end_cell;
            
            pix_sim(start_cell,:) = NaN;
            pix_sim(end_cell,:) = NaN;
        end
    end
    
    while (any(any(pix_sim > 0.2)))
        [start_cell,end_cell] = find(pix_sim == max(pix_sim(:)),1,'first');
        
        all_tracking_props{i_num}(start_cell).next_obj = end_cell;
        
        pix_sim(start_cell,:) = NaN;
        pix_sim(end_cell,:) = NaN;
        
    end    
end

%detect conflicts in pixel similarity assignment
for i_num=1:length(all_tracking_props)
    try
        predicted_next_objs = [all_tracking_props{i_num}.next_obj];
    catch %#ok<*CTCH>
        continue;
    end
    
    assert(length(predicted_next_objs) >= length(unique(predicted_next_objs)))
end

tracking_mat = convert_tracking_props_to_matrix(all_tracking_props) - 1;

output_file = fullfile(base_dir, image_dirs(1).name,filenames.tracking);

%If the tracking matrix is empty, don't output anything and don't make a
%folder for the empty matrix
if (any(size(tracking_mat) == 0))
    
else
    for track_col = 1:size(tracking_mat,2)
        track_nums = sort(unique(tracking_mat(:,1)+1));
        if (track_nums(1) == 0), track_nums = track_nums(2:end); end
        
        assert(length(track_nums) == track_nums(end));
    end
    
    if (not(exist(fileparts(output_file),'dir'))), mkdir(fileparts(output_file)); end
    
    csvwrite(fullfile(base_dir, image_dirs(1).name,filenames.tracking),tracking_mat);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tracking_matrix = convert_tracking_props_to_matrix(all_tracking_props)

cells = struct('start',{},'sequence',{});
tracking_num = 1;

[i_num,obj_num] = find_unassigned_obj(all_tracking_props);

while (i_num ~= 0 && obj_num ~= 0)
    if (length(cells) < tracking_num)
        cells(tracking_num).start = i_num;
    end
    
    cells(tracking_num).sequence = [cells(tracking_num).sequence, obj_num];
    all_tracking_props{i_num}(obj_num).assigned = 1;
    
    %pick out the next object to follow
    [i_num, obj_num] = follow_to_next_obj(all_tracking_props,i_num,obj_num);
    
    if (i_num == 0)
        [i_num, obj_num] = find_unassigned_obj(all_tracking_props);
        tracking_num = tracking_num + 1;
    end
end

tracking_matrix = zeros(length(cells),length(all_tracking_props));

for cell_num = 1:length(cells)
    col_range = cells(cell_num).start:(cells(cell_num).start + length(cells(cell_num).sequence) - 1);
    tracking_matrix(cell_num,col_range) = cells(cell_num).sequence;
end

for col_num = 1:size(tracking_matrix,2)
    this_col = tracking_matrix(:,col_num);
    this_col = this_col(this_col ~= 0);
    
    this_col = sort(unique(this_col))';
    
    %empty columns mean there weren't any adhesions in that time step, so
    %check for that in the following assert first, then if there were
    %cells, make sure all were accounted for
    assert(isempty(this_col) || all(this_col == 1:max(this_col)))
    
    assert((isempty(all_tracking_props{col_num}) && isempty(this_col)) || ...
        (length(all_tracking_props{col_num}) == max(this_col)));
end

function [i_num,obj_num] = find_unassigned_obj(all_tracking_props)

for i_num=1:length(all_tracking_props)
    for obj_num=1:size(all_tracking_props{i_num},1)
        try
            if (all_tracking_props{i_num}(obj_num).assigned)
                continue;
            else
                return;
            end
        catch
            continue;
        end
    end
end

i_num = 0;
obj_num = 0;

return;

function [i_num,obj_num] = follow_to_next_obj(all_tracking_props,i_num,obj_num)

try %#ok<*TRYNC>
    if (size(all_tracking_props{i_num}(obj_num).next_obj,2) == 0)
        i_num = 0;
        obj_num = 0;
        return;
    else
        obj_num = all_tracking_props{i_num}(obj_num).next_obj;
        i_num = i_num + 1;
        return;
    end
end

i_num = 0;
obj_num = 0;

return;