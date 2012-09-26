function perim_thick = thicken_perimeter(perims,objects,varargin)
%THICKEN_PERIMETER    Thicken a provided set of perimeters
%
%   thick_perims = thicken_perimeter(thin_perims,objects) thickens the
%   perimeters in thin_perims to 3 pixels, where the thickening is only
%   allowed within the matching object provided in objects
%
%   Options:
%       -3rd argument: thickeness of perimeter, defaults to 3

if (islogical(perims)), perims = double(perims); end
if (islogical(objects)), 
    objects = double(objects);
else
    objects = ones(size(perims));
end

thickness = 3;
if (size(varargin,2) > 0)
    thickness = varargin{1};
end

perim_thick = perims;
three_disk = strel('disk',thickness,0);
for perim_num = 1:max(perims(:))
    this_cell_perim = perims == perim_num;
    cell_nums = unique(objects(this_cell_perim));
    
    this_cell = zeros(size(objects,1),size(objects,2));
    for cell_num = cell_nums'
        this_cell(objects == cell_num) = 1;
    end
     
    this_cell_perim = imdilate(this_cell_perim,three_disk);
    this_cell_perim = this_cell_perim & this_cell;
    
    perim_thick(this_cell_perim) = perim_num;
end
