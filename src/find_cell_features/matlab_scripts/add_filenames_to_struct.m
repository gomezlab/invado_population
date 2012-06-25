function filenames = add_filenames_to_struct(filenames)
%ADD_FILENAMES_TO_STRUCT    Takes in a struct, to which the filenames for
%                          various files made in the pipeline are added

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Main Program
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filenames.puncta_labeled = 'puncta_labeled.png';
filenames.puncta = 'puncta.png';
filenames.puncta_median = '../../../all_field_resources/puncta_median.png';
filenames.puncta_range = '../../../all_field_resources/puncta_image_range.csv';

filenames.cell_mask = 'cell_mask.png';
filenames.highlighted_cell_mask = 'highlighted_mask.png';
filenames.labeled_cell_mask = 'cell_mask_labeled.png';
filenames.labeled_cell_mask_perim = 'cell_mask_perim_labeled.png';

filenames.gel = 'gel.png';
filenames.gel_median = '../../../all_field_resources/gel_median.png';
filenames.gel_range = '../../../all_field_resources/gel_image_range.csv';

filenames.intensity_correction = 'intensity_correction.csv';

filenames.no_cells = '../../cell_props/no_cell_regions.png';

filenames.invader_vis = 'invader_vis.png';

filenames.tracking_raw = '../../cell_props/raw_tracking_data.mat';
filenames.tracking = '../../tracking_matrices/tracking_seq.csv';
filenames.tracking_vis = 'tracking.png';

filenames.cell_props = '../../cell_props/cell_props.mat';
filenames.active_degrade = '../../cell_props/active_degrade.csv';

filenames.lineage_dir = '../../cell_props/lin_time_series/';

filenames.single_cell_dir = '../../cell_props/single_cells/';
filenames.final_gel_diffs = '../../cell_props/final_gel_diffs.csv';
filenames.corrected_final_gel_diffs = '../../cell_props/corrected_final_gel_diffs.csv';