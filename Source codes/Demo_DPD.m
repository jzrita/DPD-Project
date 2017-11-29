clear
clc
load('.\colormap.mat');
%% ================== generate seed  =======================
i=1;
data_dir='.\SR4000_test_images';
im_base_name = sprintf( 'Scene%d', i );
im_file = sprintf( '%s\\%s.bmp', data_dir, im_base_name );
im_d = double(imread(im_file));
seeds_data_in.initial_seeds_shape = ones(4);
seeds_data_in.seed_shape = 'square';
seeds_data_in.order_of_generating_seeds = 'RasterScan';
seeds_data_in.im_d = im_d;
seeds_data_in.places_to_cehck_for_initial_seeding = im_d;
seeds_data_in.min_allowed_points_for_initial_seeding = sum(sum(seeds_data_in.initial_seeds_shape));
seeds_data = get_initial_seeds(seeds_data_in,im_base_name);
%% ================== plane growing  =======================
fields_data_in.seeds_data{1}.im_d = [];
fields_data_in.im_file = im_file;
fields_data_in.im_d = im_d;
fields_data_in.places_need_to_be_cehcked_for_seeding = ones(size(fields_data_in.im_d));
fields_data_in.order_of_using_seeds = 'ascendorder';
fields_data_in.basic_data.lag_of_update_THR = 1;
fields_data_in.basic_data.maximum_allowed_dissimilarity = 3;
fields_data_in.seeds_data{1} = seeds_data;
weight=0.009;

[fields_data, ~] = region_growing(fields_data_in,weight);
figure;image(uint8(fields_data.fields));colormap(settelments_colormap);axis off;
%% ================== over-grwoing correction  =======================
epoch_idx = 1;
basic_penetrating_element_shape = 'line';
show.show_im = 1;
show.fig_idx = 2;
show.settelments_colormap = settelments_colormap;
[fields_data_out] = over_growing_correction(fields_data, basic_penetrating_element_shape, epoch_idx, show);
figure;image(uint8(fields_data_out.fields));colormap(settelments_colormap);

%% ================== under-grwowing correction  =======================
fields_data_out1 =  under_growing_correction(fields_data_out);
figure;image(uint8(fields_data_out1.parallel_surface_detection.field_index));colormap(settelments_colormap); axis off;

%% ================== calculate ROC  =======================

if i==1
%%%%%%%%%%%%%%%%%%% Scene 1 %%%%%%%%%%%%%%%%%%%%%%%%%
im_GT_name='Scene1_GT.bmp';
Plane_index_gt=[129 225 213 86 226 88 116 35];
elseif i==2
%%%%%%%%%%%%%%%%%%% Scene 2 %%%%%%%%%%%%%%%%%%%%%%%%%
 im_GT_name='Scene2_GT_new.bmp';
 Plane_index_gt=[239 55 232 213 62];
elseif i==3
%%%%%%%%%%%%%%%%%%% Scene 3 %%%%%%%%%%%%%%%%%%%%%%%%%
 im_GT_name='Scene3_GT.bmp';
 Plane_index_gt=[83 155 123 218 180 18 225];
elseif i==4
 %%%%%%%%%%%%%%%%%%% Scene 4 %%%%%%%%%%%%%%%%%%%%%%%%%
 im_GT_name='Scene4_GT.bmp';
 Plane_index_gt=[123 155 137 200 201 217];
elseif i==5
 %%%%%%%%%%%%%%%%%%% Scene 5 %%%%%%%%%%%%%%%%%%%%%%%%%
 im_GT_name='Scene5_GT_new1.bmp';
 Plane_index_gt=[239 210 251 135 103 7 113 213 62 102];
end

check_threshold=0.8;
im_gt=imread(im_GT_name);
if size(im_gt,3)>1
im_gt=rgb2gray(im_gt);
end
im_test=uint8(fields_data_out1.parallel_surface_detection.field_index);
%im_test=uint8(fields_data_out.fields);
%im_test=uint8(fields_data.fields);
[~,~,~,Sensitivity,correct_detected_plane,Specificity]=calculateROC(im_gt,im_test,check_threshold,Plane_index_gt);
Sensitivity_mean=mean(Sensitivity)
Specificity_mean=mean(Specificity)
CDR=(sum(correct_detected_plane)/length(Plane_index_gt))*100


