function diffs = collect_local_diff_properties(data_struct,this_ad,varargin)
% COLLECT_LOCAL_DIFF_PROPERTIES    find the local difference properties of
%                                  a provided set of images and a given
%                                  puncta

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Setup variables and parse command line
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
i_p = inputParser;
i_p.FunctionName = 'COLLECT_LOCAL_DIFF_PROPERTIES';

i_p.addRequired('data_struct',@isstruct);
i_p.addRequired('this_ad',@islogical);

i_p.addParamValue('background_border_size',5,@(x)isnumeric(x));
i_p.addOptional('debug',0,@(x)x == 1 || x == 0);

i_p.parse(data_struct,this_ad,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

diffs = struct();

background_region = logical(imdilate(this_ad,strel('disk',i_p.Results.background_border_size,0)));
%we don't want to include any areas that have been identified as other
%objects
background_region = and(background_region,not(data_struct.adhesions));
%also exclude areas outside the registered image
background_region = logical(background_region .* data_struct.binary_shift);

if (sum(sum(background_region)) == 0)
    diffs.Local_gel_diff = NaN;
    diffs.Local_gel_diff_corr = NaN;
else
    diffs.Local_gel_diff = mean(data_struct.gel_image(this_ad)) - mean(data_struct.gel_image(background_region));
    diffs.Local_gel_diff_corr = diffs.Local_gel_diff*data_struct.intensity_correction;
end
diffs.Global_gel_diff = mean(data_struct.gel_image(this_ad)) - mean(data_struct.gel_image(data_struct.binary_shift));

large_background_region = logical(imdilate(this_ad,strel('disk',10,0)));
%we don't want to include any areas that have been identified as other
%objects
large_background_region = and(large_background_region,not(data_struct.adhesions));
%also exclude areas outside the registered image
large_background_region = logical(large_background_region .* data_struct.binary_shift);

if (sum(sum(large_background_region)) == 0)
    diffs.Large_local_gel_diff = NaN;
    diffs.Large_local_gel_diff_corr = NaN;
else
    diffs.Large_local_gel_diff = mean(data_struct.gel_image(this_ad)) - mean(data_struct.gel_image(large_background_region));
    diffs.Large_local_gel_diff_corr = diffs.Large_local_gel_diff*data_struct.intensity_correction;
end