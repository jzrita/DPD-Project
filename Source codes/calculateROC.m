

function [count_TP,count_FN,count_FP,Sensitivity,correct_detected_plane,Specificity]=calculateROC_NYU(im_gt,im_test,check_threshold,Plane_index_gt)
%function [count_TP,count_FN,count_FP,Sensitivity,correct_detected_plane,Specificity]=calculateROC(im_GT_name,im_test_name,check_threshold,Plane_index_gt)
% im_gt=imread(im_GT_name);
% %im_gt=rgb2gray(im_gt);
% figure;imshow(im_gt);
%
% im=imread(im_test_name);
% im=rgb2gray(im);
% im_test=imresize(im,size(im_gt),'nearest');
% figure;imshow(im_test);

correct_detected_plane=zeros(1,length(Plane_index_gt));
count_TP=zeros(1,length(Plane_index_gt));
count_TN=zeros(1,length(Plane_index_gt));
count_FN=zeros(1,length(Plane_index_gt));
count_FP=zeros(1,length(Plane_index_gt));
Sensitivity=zeros(1,length(Plane_index_gt));
Specificity=zeros(1,length(Plane_index_gt));
plane_index_have_been_checked=zeros(1,length(Plane_index_gt));
count=1;
for i=1:length(Plane_index_gt)
    [r_gt,c_gt]=find(im_gt==Plane_index_gt(i));
    index=sub2ind(size(im_gt),r_gt,c_gt);
    plane_point=im_test(index);
    plane_point_indensity=unique(plane_point,'rows'); %find all the plane indexs in this area, delete the repeat one
    plane_point_indensity(plane_point_indensity==0) = [];%delete holes
    points_intensity_number=zeros(1,length(plane_point_indensity));
    for j=1:length(plane_point_indensity)
        points_intensity_number(j)=sum(plane_point==plane_point_indensity(j)); %calculate the points number for each plane index
    end
    [~,plane_index_in_order]=sort(points_intensity_number,'descend');
    %idx_max=find(points_intensity_number==max(points_intensity_number)); %find the max one as the main plane is this area
    idx_max=plane_index_in_order(1);
    if isempty(find(plane_index_have_been_checked==plane_point_indensity(idx_max),1))
        plane_index_have_been_checked(count)=plane_point_indensity(idx_max);
        count=count+1;
        count_TP(i)=sum(plane_point==plane_point_indensity(idx_max)); % the correct points in this area
        count_FN(i)=length(plane_point)-count_TP(i); % the points belong to other planes in this area
        Sensitivity(i)=(count_TP(i)/(count_TP(i) + count_FN(i)))*100;
        if count_TP(i)>=check_threshold*length(plane_point)
            correct_detected_plane(i)=1;
        end
        im_test(index)=0;
        [r_edge,c_edge]=find(im_gt==255);
        index_edge=sub2ind(size(im_gt),r_edge,c_edge);
        im_test(index_edge)=0;
        [r_max,~]=find(im_test==plane_point_indensity(idx_max));
        if ~isempty(r_max)
            count_FP(i)=length(r_max); %in other areas find the points belong to this plane
        end
        count_TN(i)=size(im_gt,1)*size(im_gt,2)-count_FP(i)-length(plane_point); %other plane points in other areas
        Specificity(i)=count_TN(i)/(count_TN(i)+count_FP(i))*100;
    else
        % [~, index_tem]=setxor(plane_point_indensity,plane_point_indensity(idx_max));
        % plane_points_intensity_new=plane_point_indensity(sort(index_tem)); %delete the already checked plane index
        % idx_max_second=find(plane_points_intensity_new==max(plane_points_intensity_new)); %find the next max one as the main plane is this area
       if length(plane_index_in_order)>2
        idx_max_second= plane_index_in_order(2);
            [r_tem,c_tem]=find(im_test==plane_point_indensity(idx_max_second)); %find all the points that belong to the same plane
            index_tem1=sub2ind(size(im_test),r_tem,c_tem);
            size_of_second_max_plane=sum(plane_point==plane_point_indensity(idx_max_second));
            if size_of_second_max_plane >10 && size_of_second_max_plane/length(index_tem1)>0.5 %to check whether most part of the plane is in this area
                count_TP(i)=sum(plane_point==plane_point_indensity(idx_max_second)); % the correct points in this area
                count_FN(i)=length(plane_point)-count_TP(i); % the points belong to other planes in this area
                Sensitivity(i)=(count_TP(i)/(count_TP(i) + count_FN(i)))*100;
                if count_TP(i)>=check_threshold*length(plane_point)
                    correct_detected_plane(i)=1;
                end
                im_test(index)=0;
                [r_edge,c_edge]=find(im_gt==255);
                index_edge=sub2ind(size(im_gt),r_edge,c_edge);
                im_test(index_edge)=0;
                [r_second,~]=find(im_test==plane_point_indensity(idx_max_second));
                if ~isempty(r_second)
                    count_FP(i)=length(r_second); %in other areas find the points belong to this plane
                end
                count_TN(i)=size(im_gt,1)*size(im_gt,2)-count_FP(i)-length(plane_point);
                Specificity(i)=count_TN(i)/(count_TN(i)+count_FP(i))*100;
            else
                idx_max_third= plane_index_in_order(3);
                if ~isempty(idx_max_third)
                    [r_tem,c_tem]=find(im_test==plane_point_indensity(idx_max_third)); %find all the points that belong to the same plane
                    index_tem1=sub2ind(size(im_test),r_tem,c_tem);
                    size_of_third_max_plane=sum(plane_point==plane_point_indensity(idx_max_third));
                    if size_of_third_max_plane >10 && size_of_third_max_plane/length(index_tem1)>0.5 %to check whether most part of the plane is in this area
                        count_TP(i)=sum(plane_point==plane_point_indensity(idx_max_third)); % the correct points in this area
                        count_FN(i)=length(plane_point)-count_TP(i); % the points belong to other planes in this area
                        Sensitivity(i)=(count_TP(i)/(count_TP(i) + count_FN(i)))*100;
                        if count_TP(i)>=check_threshold*length(plane_point)
                            correct_detected_plane(i)=1;
                        end
                        im_test(index)=0;
                        [r_edge,c_edge]=find(im_gt==255);
                        index_edge=sub2ind(size(im_gt),r_edge,c_edge);
                        im_test(index_edge)=0;
                        [r,~]=find(im_test==plane_point_indensity(idx_max_third));
                        if ~isempty(r)
                            count_FP(i)=length(r); %in other areas find the points belong to this plane
                        end
                        count_TN(i)=size(im_gt,1)*size(im_gt,2)-count_FP(i)-length(plane_point);
                        Specificity(i)=count_TN(i)/(count_TN(i)+count_FP(i))*100;
                    else
                        Sensitivity(i)=0;
                        Specificity(i)=0;
                    end
                else
                    Sensitivity(i)=0;
                    Specificity(i)=0;
                end
            end
        else
            Sensitivity(i)=0;
            Specificity(i)=0;
        end
    end
end

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% Scene 1 %%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc
% im_GT_name='Scene1_GT.bmp';
% Plane_index_gt=[129 225 213 86 226 88 116 35];
% im_test_name='CC-RANSAC_table_chair.png';%'RPCA_table_chair.png';%'CORG_table_chair.png';%
%%%%%%%%%%%%%%%%%%%% Scene 2 %%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc
% im_GT_name='Scene2_GT_new.bmp';
% Plane_index_gt=[239 55 232 213 62];
% im_test_name='RPCA_table_seat.png';%'CORG_table_seat.png';%'CC-RANSAC_table_seat.png';%
%%%%%%%%%%%%%%%%%%% Scene 3 %%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc
% im_GT_name='Scene3_GT.bmp';
% Plane_index_gt=[83 155 123 218 180 18 225];
% im_test_name='CORG_cabinet.png';%'RPCA_cabinet.png';%'CC-RANSAC_cabinet.png';%
% %%%%%%%%%%%%%%%%%%% Scene 4 %%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc
% im_GT_name='Scene4_GT.bmp';
% Plane_index_gt=[123 155 137 200 201 217];
% im_test_name='CC-RANSAC_stair1.png';%'CORG_stair1.png';%'RPCA_stair1.png';%
% %%%%%%%%%%%%%%%%%%% Scene 5 %%%%%%%%%%%%%%%%%%%%%%%%%
% clear;
% clc
% im_GT_name='Scene5_GT_new1.bmp';
% Plane_index_gt=[239 210 251 135 103 7 113 102];%213 62 102];
% im_test_name='RPCA_stair2.png';%'CORG_stair2.png';%'CC-RANSAC_stair2.png';%
%
% check_threshold=0.8;
% [count_TP,count_FN,count_FP,Sensitivity,correct_detected_plane,Specificity]=calculateROC(im_GT_name,im_test_name,check_threshold,Plane_index_gt);
% mean(Sensitivity)
% mean(Specificity)
% (sum(correct_detected_plane)/length(Plane_index_gt))*100
%(Precise_rate>threshold_correct_detected) = correct_detected_plane;
% Recall_rate indicates the over-growing rate, the larger the better(less over growing);
% FN indicates the under-growing(over-segment);
% FP indicates the over-growing;

% load('D:\Year1\Project_0_Resubmission_of_DDPSD\Test_Results\Indoor_S1-1_square_016_ascendorder_after_PP.mat');
% load('D:\Year1\Project_0_Resubmission_of_DDPSD\Test_Results\Indoor_S1-7_square_016_ascendorder_after_PP.mat');
% load('D:\Year1\Project_0_Resubmission_of_DDPSD\Test_Results\Indoor_S2-1_square_016_ascendorder_after_PP.mat');
% load('D:\Year1\Project_0_Resubmission_of_DDPSD\Test_Results\stair1_square_016_ascendorder_after_PP.mat');
%  load('D:\Year1\Project_0_Resubmission_of_DDPSD\Test_Results\stair2_square_016_ascendorder_after_PP.mat');
%  im_test=uint8(fields_data_out.parallel_surface_detection.field_index);