function normed_image = normalize_image(image,varargin)

if (length(varargin) >= 1)
    i_range = varargin{1};
else
    i_range = [min(image(:)),max(image(:))];
end

normed_image = (image - i_range(1))/range(i_range);