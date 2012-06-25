function build_all_montage_file_sets(exp_dir,varargin)
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

build_raw_montage_file_set(exp_dir,'gel.png')
build_raw_montage_file_set(exp_dir,'focal_image.png')
build_raw_montage_file_set(exp_dir,'highlighted_mask.png')
build_raw_montage_file_set(exp_dir,'cell_overlaps.png')
build_raw_montage_file_set(exp_dir,'uncorrected_focal_image.png')
build_raw_montage_file_set(exp_dir,'gel_invader.png')
build_raw_montage_file_set(exp_dir,'tracking.png')
