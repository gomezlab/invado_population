function check_degradation_gui(exp_dir,varargin)
tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;

i_p.addRequired('exp_dir',@(x)exist(x,'dir') == 7);

i_p.addParamValue('debug',0,@(x)x == 1 || x == 0);

i_p.parse(exp_dir,varargin{:});

addpath(genpath('..'));

global file_data;
file_data = struct();
file_data.exp_dir = exp_dir;
file_data.filenames = add_filenames_to_struct(struct());

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global i_num;
i_num = 1;

global field_num;
field_num = 1;

file_data.fields = dir(file_data.exp_dir);
file_data.fields = filter_to_time_series(file_data.fields);
load_field_dir();

global all_cmaps;
all_cmaps = cell(length(file_data.image_dirs),length(file_data.fields));

global images;
images = struct();
load_new_images();

close all;

figure

global ui_controls;
ui_controls = struct();

ui_controls.status = uicontrol('Style','text','Units','normalized','Position',[0.01,0.80,0.6,0.15], ... 
    'FontSize',16);
update_status();

uicontrol('Style','pushbutton','String','Prev Field', ... 
    'Units','normalized','Position',[0.7,0.8,0.14,0.15], ...
    'Callback',@jump_to_prev_field);
uicontrol('Style','pushbutton','String','Next Field', ... 
    'Units','normalized','Position',[0.85,0.8,0.14,0.15], ...
    'Callback',@jump_to_next_field);

ax(1) = subplot(1,2,1);
ui_controls.high_hnd = image(images.gel_high,'ButtonDownFcn',@get_clicks);

ax(2) = subplot(122);
ui_controls.no_high_hnd = image(images.gel);

linkaxes(ax,'xy')
axis(ax,'image')

uicontrol('Style','pushbutton','String','Image 1', ... 
    'Units','normalized','Position',[0.01,0.01,0.10,0.15], ...
    'Callback',@jump_to_beginning);

uicontrol('Style','pushbutton','String','Previous Image', ... 
    'Units','normalized','Position',[0.11,0.01,0.5-0.11,0.15], ...
    'Callback',@jump_to_previous);
uicontrol('Style','pushbutton','String','Next Image',...
    'Units','normalized','Position',[0.5,0.01,0.39,0.15], ...
    'Callback',@jump_to_next);

uicontrol('Style','pushbutton','String','Play Movie', ... 
    'Units','normalized','Position',[0.89,0.01,0.10,0.15], ...
    'Callback',@play_movie);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cell_num = search_cell_hits(x,y,props)
x = round(x);
y = round(y);

cell_num = 0;
for i=1:length(props)
    cell_pixels = props(i).PixelList;
    
    x_hits = cell_pixels(:,1) == x;
    y_hits = cell_pixels(:,2) == y;
    
    if (any(x_hits & y_hits))
        cell_num = i;
    end
end

function get_clicks(src,eventdata) %#ok<*INUSD>
x = 1;
y = 1;

global images;
global all_cmaps;
global i_num;
global field_num;
global ui_controls;

while(x > 0 && y > 0) 

    [x,y,button]=ginput(1);
    
    if (button ~= 1)
        break;
    end
    
    cell_num = search_cell_hits(x,y,images.props);
    if (cell_num > 0)
        if (all(all_cmaps{i_num,field_num}(cell_num,:) == [0,1,0]))
            all_cmaps{i_num,field_num}(cell_num,:) = [1,0,0];
        else
            all_cmaps{i_num,field_num}(cell_num,:) = [0,1,0];
        end
        
        new_highlight = create_highlighted_image(images.gel,images.cell_perim,'color_map',all_cmaps{i_num,field_num});
        set(ui_controls.high_hnd,'CData',new_highlight);
    end
end

function load_new_images()

global images;
global file_data;
global all_cmaps;
global i_num;
global field_num;

images.cell_label = imread(fullfile(file_data.field_dir,file_data.image_dirs(i_num).name,file_data.filenames.labeled_cell_mask_filename));
images.cell_perim = imread(fullfile(file_data.field_dir,file_data.image_dirs(i_num).name,file_data.filenames.labeled_cell_mask_perim));

images.gel = imread(fullfile(file_data.field_dir,file_data.image_dirs(i_num).name,file_data.filenames.gel_filename));
images.gel_range = csvread(fullfile(file_data.field_dir,file_data.image_dirs(i_num).name,file_data.filenames.gel_range_file));

images.gel = (double(images.gel) - images.gel_range(1))/range(images.gel_range);
images.gel = cat(3,images.gel,images.gel,images.gel);

images.props = regionprops(images.cell_label,'PixelList');

if (isempty(all_cmaps{i_num,field_num}))
    all_cmaps{i_num,field_num} = zeros(max(images.cell_label(:)),3);
    all_cmaps{i_num,field_num}(:,1) = 1;
end

images.gel_high = create_highlighted_image(images.gel,images.cell_perim,'color_map',all_cmaps{i_num,field_num});

function jump_to_next(src,eventdata)

global i_num;
global file_data;
global images;

global ui_controls;

if ((i_num + 1) <= length(file_data.image_dirs))
    i_num = i_num + 1;
else
    i_num = 1;
end

load_new_images();
set(ui_controls.high_hnd,'CData',images.gel_high);
set(ui_controls.no_high_hnd,'CData',images.gel);
update_status();

function jump_to_next_field(src,eventdata)

global i_num;
global field_num;
global file_data;
global images;

global ui_controls;

if ((field_num + 1) <= length(file_data.fields))
    field_num = field_num + 1;
else
    field_num = 1;
end

i_num = 1;

load_field_dir();
load_new_images();
set(ui_controls.high_hnd,'CData',images.gel_high);
set(ui_controls.no_high_hnd,'CData',images.gel);
update_status();

function load_field_dir

global field_num;
global file_data;

file_data.field_dir = fullfile(file_data.exp_dir,file_data.fields(field_num).name,'individual_pictures');
file_data.image_dirs = dir(file_data.field_dir);

assert(strcmp(file_data.image_dirs(1).name, '.'), 'Error: expected "." to be first string in the dir command')
assert(strcmp(file_data.image_dirs(2).name, '..'), 'Error: expected ".." to be second string in the dir command')
assert(str2num(file_data.image_dirs(3).name) == 1, 'Error: expected the third string to be image set one') %#ok<ST2NM>

file_data.image_dirs = file_data.image_dirs(3:end);

function jump_to_previous(src,eventdata)

global i_num;
global file_data;
global images;

global ui_controls;

if ((i_num - 1) > 0)
    i_num = i_num - 1;
else
    i_num = length(file_data.image_dirs);
end

load_new_images();
set(ui_controls.high_hnd,'CData',images.gel_high);
set(ui_controls.no_high_hnd,'CData',images.gel);
update_status();

function jump_to_prev_field(src,eventdata)

global i_num;
global field_num;
global file_data;
global images;

global ui_controls;

if ((field_num - 1) > 0)
    field_num = field_num - 1;
else
    field_num = length(file_data.fields);
end

i_num = 1;

load_field_dir();
load_new_images();
set(ui_controls.high_hnd,'CData',images.gel_high);
set(ui_controls.no_high_hnd,'CData',images.gel);
update_status();

function play_movie(src,eventdata)

global i_num;
global file_data;

for i=i_num:(length(file_data.image_dirs)-1)
    
    jump_to_next(src,eventdata);
	drawnow;
end

function jump_to_beginning(src,eventdata)

global i_num;
global images;

global ui_controls;

i_num = 1;

load_new_images();
set(ui_controls.high_hnd,'CData',images.gel_high);
set(ui_controls.no_high_hnd,'CData',images.gel);
update_status();

function update_status()

global i_num;
global field_num;
global file_data;

global ui_controls;

status_str = {['Field: ',file_data.fields(field_num).name], ...
    ['Image Number: ', num2str(i_num), '/', num2str(length(file_data.image_dirs))]};
    

set(ui_controls.status,'String',status_str);