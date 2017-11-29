

function fields_data= under_grwowing_correction(fields_data)
global number_of_plane
global r
global c
global NumberofEmptyField
number_of_plane=length(fields_data.field_);
[r,c]=size(fields_data.fields);
NumberofEmptyField=0;
for i=1:number_of_plane
    if isempty(fields_data.field_{i})
        NumberofEmptyField=NumberofEmptyField+1;
        continue
    else
        size_of_surface(i)=length(fields_data.field_{i}.points_of_filed);
    end
end
%******************************************************************************************
[~, order_of_process] = sort(size_of_surface, 'descend');

% Extract the normal vector for each surface from largest surface to smallest one and Calculate all the angles of vectors of one surface
[field_index,fields_data.field_, max_ang_of_plane_matrix]=get_normal_vect_and_angles(fields_data.field_,order_of_process);

colored_field_index = zeros(r, c, 3);
if max(max(field_index)) < 0.95 * 2^24
    tmp_im = 0.95 * 2^24 * ( field_index / max(max(field_index)) );
    tmp_im = uint64(round(tmp_im -10));
    colored_field_index(:, :, 1) = bitand(tmp_im,  2^8-1); %bit-wise AND; bit by bit AND
    colored_field_index(:, :, 2) = bitand(bitshift(tmp_im, -8),  2^8-1);%(jz) bitshift(tmp_im, -8), the bits in tmo_im shift to right by 8 bits
    colored_field_index(:, :, 3) = bitand(bitshift(tmp_im, -16),  2^8-1);
end
%  Calculate the angle differences between two different surface
ang_of_two_planes_matrix=calculate_angle_differences(fields_data.field_,order_of_process);

% Detect the parallel surfaces
[parallel_plane,parallel_surface_matrix_down]=get_parallel_plane(fields_data.field_,ang_of_two_planes_matrix,max_ang_of_plane_matrix,order_of_process);

% Detect the aligned surface
[aligned_plane,times_of_alignment,im_deviation_of_points,fields_data.field_]=detect_aligned_plane(fields_data.field_,order_of_process, parallel_plane);

% merged surface
[field_index,merged_plane,times_of_merge,field_delete_merged_plane]=merge_surface(fields_data.field_, aligned_plane,field_index);

count=1;
need_to_delete=[];
for i=1:number_of_plane
    if isempty(field_delete_merged_plane{i})
        need_to_delete(count)=i;
        continue;
    else
        if isempty(field_delete_merged_plane{i}.points_of_filed)
            need_to_delete(count)=i;
            count=count+1;
        end
    end
end
if  isempty(need_to_delete)~=1
    field_delete_merged_plane(need_to_delete) = [];
end

merged_colored_fields = zeros(r, c, 3);
if max(max(field_index)) < 0.95 * 2^24
    tmp_im = 0.95 * 2^24 * ( field_index / max(max(field_index)) );
    tmp_im = uint64(round(tmp_im -10));
    merged_colored_fields(:, :, 1) = bitand(tmp_im,  2^8-1); %bit-wise AND; bit by bit AND
    merged_colored_fields(:, :, 2) = bitand(bitshift(tmp_im, -8),  2^8-1);%(jz) bitshift(tmp_im, -8), the bits in tmo_im shift to right by 8 bits
    merged_colored_fields(:, :, 3) = bitand(bitshift(tmp_im, -16),  2^8-1);
end

%%%%%%%%%%%% Output %%%%%%%%%%%%%%%%
fields_data.parallel_surface_detection.parallel_surface_matrix_down=parallel_surface_matrix_down;
fields_data.parallel_surface_detection.field_index=field_index;
fields_data.parallel_surface_detection.field_delete_merged_plane=field_delete_merged_plane;
fields_data.parallel_surface_detection.colored_field_index=colored_field_index;
fields_data.parallel_surface_detection.parallel_plane=parallel_plane;
fields_data.parallel_surface_detection.im_deviation_of_points=im_deviation_of_points;
if times_of_merge>1
    fields_data.parallel_surface_detection.merged_colored_fields=merged_colored_fields;
    fields_data.parallel_surface_detection.merged_plane=merged_plane;
end
if times_of_alignment>1
    fields_data.parallel_surface_detection.aligned_plane=aligned_plane;
end

return


function [field_index,field,max_ang_of_plane_matrix]=get_normal_vect_and_angles(field,order_of_process)
global number_of_plane
global r
global c
global NumberofEmptyField
max_ang_of_plane=zeros(1,number_of_plane-NumberofEmptyField);
r_points=cell(1,number_of_plane-NumberofEmptyField);
c_points=cell(1,number_of_plane-NumberofEmptyField);
for i=1:(number_of_plane-NumberofEmptyField)
    r_points{i}(:,1)=field{order_of_process(i)}.points_of_filed(:,1);
    c_points{i}(:,1)=field{order_of_process(i)}.points_of_filed(:,2);
end
%******************************************************************************************
field_index=zeros(r,c);

for idx_plane=1:(number_of_plane-NumberofEmptyField)
    field_index(sub2ind([r,c], r_points{idx_plane}, c_points{idx_plane}))=order_of_process(idx_plane);
    % Extract the normal vector for each surface from largest surface to smallest one
    for idx_equation=1:length(field{order_of_process(idx_plane)}.equations_of_plan_)
        field{order_of_process(idx_plane)}.Normal_vector(idx_equation,:)=field{order_of_process(idx_plane)}.equations_of_plan_{idx_equation};
    end
    % Calculate all the angles of vectors of one surface
    field{order_of_process(idx_plane)}.Normal_vector_final=field{order_of_process(idx_plane)}.Normal_vector(end,:);
    if length(field{order_of_process(idx_plane)}.equations_of_plan_)>1
        for idx_equation1=1:length(field{order_of_process(idx_plane)}.equations_of_plan_)-1
            field{order_of_process(idx_plane)}.angle_of_plane(idx_equation1)=calculate_angle(field{order_of_process(idx_plane)}.Normal_vector_final(1,1:3), field{order_of_process(idx_plane)}.Normal_vector(idx_equation1,1:3));
        end
        max_ang_of_plane(idx_plane)=max(field{order_of_process(idx_plane)}.angle_of_plane);
    else  max_ang_of_plane(idx_plane)=0;
    end
end
max_ang_of_plane_matrix=repmat(max_ang_of_plane,number_of_plane-NumberofEmptyField,1);
return%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function ang_of_two_planes_matrix=calculate_angle_differences(field,order_of_process)
global number_of_plane
global NumberofEmptyField
ang_of_two_planes_matrix=zeros(number_of_plane-NumberofEmptyField);
for ii=1:(number_of_plane-NumberofEmptyField)-1
    for jj=ii+1:(number_of_plane-NumberofEmptyField)
        ang_of_two_planes_matrix(jj,ii)=calculate_angle(field{order_of_process(ii)}.Normal_vector_final(1,1:3),field{order_of_process(jj)}.Normal_vector_final(1,1:3));
    end
end
return%||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function ang_of_plane=calculate_angle(vect1,vect2)
A1=vect1(1);
B1=vect1(2);
C1=vect1(3);
A2=vect2(1);
B2=vect2(2);
C2=vect2(3);
ang_of_plane=(180/pi)*acos(sum( [A1, B1,C1].*[A2, B2,C2])/((sqrt(A1^2+B1^2+C1^2))*(sqrt(A2^2+B2^2+C2^2))));
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [parallel_plane,parallel_surface_matrix_down]=get_parallel_plane(field,ang_of_two_planes_matrix,max_ang_of_plane_matrix,order_of_process)
parallel_surface_matrix=ang_of_two_planes_matrix < max_ang_of_plane_matrix;
parallel_surface_matrix_down=tril(parallel_surface_matrix); %set the upper triangle matrix to zero;
parallel_surface_matrix_down=parallel_surface_matrix_down-diag(diag(parallel_surface_matrix_down)); %set the main diagornal element to zero
[parallel_r_idx,parallel_c_idx]=find(parallel_surface_matrix_down);
parallel_plane(:,1)=order_of_process(parallel_r_idx);
parallel_plane(:,2)=order_of_process(parallel_c_idx);
original_length=length(parallel_plane);
for ii=1:original_length-1
    %(this will caused some repeat like 1//2 and 2//1)
    for jj=ii+1:original_length
        if ~isempty(intersect(parallel_plane(ii,:),parallel_plane(jj,:)))
            temp_plane=setxor(parallel_plane(ii,:),parallel_plane(jj,:));
            if length(field{temp_plane(1)}.points_of_filed)>length(field{temp_plane(2)}.points_of_filed)
                parallel_plane(end+1,:)= [temp_plane(2),temp_plane(1)];
            else
                parallel_plane(end+1,:)= [temp_plane(1),temp_plane(2)];
            end
        end
    end
end
parallel_plane=unique(parallel_plane,'rows');  %delete the repeated rows

return%||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [aligned_plane,times_of_alignment,im_deviation_of_points,field]=detect_aligned_plane(field,order_of_process, parallel_plane)
global number_of_plane
global r
global c
global NumberofEmptyField
times_of_alignment=1;
im_deviation_of_points=zeros(r,c);
threshold_align_plane=0.75;

for idx_plane =1 :(number_of_plane-NumberofEmptyField)
    [field{order_of_process(idx_plane)}.deviation_of_points,~,field{order_of_process(idx_plane)}.deviation_of_each_point_vct]=measure_points_distances_from_plane(field{order_of_process(idx_plane)}.points_of_filed,field{order_of_process(idx_plane)}.Normal_vector_final);
    im_deviation_of_points(sub2ind([r,c], field{order_of_process(idx_plane)}.points_of_filed(:,1), field{order_of_process(idx_plane)}.points_of_filed(:,2))) =field{order_of_process(idx_plane)}.deviation_of_points;
end

for i=1:length(parallel_plane)
    [~,small2big_sum_of_relative_distance,small2big_relative_distance_vct] =measure_points_distances_from_plane(field{parallel_plane(i,1)}.points_of_filed, field{parallel_plane(i,2)}.Normal_vector_final);
    big2big_relative_distance_vct=field{parallel_plane(i,2)}.deviation_of_each_point_vct;
    max_small2big_relative_distance=max(small2big_relative_distance_vct);
    min_small2big_relative_distance=min(small2big_relative_distance_vct);
    max_big2big_relative_distance=max(big2big_relative_distance_vct);
    min_big2big_relative_distance=min(big2big_relative_distance_vct);
    if min_small2big_relative_distance < max_big2big_relative_distance
        xbins=linspace(min_small2big_relative_distance,max_big2big_relative_distance,100);
    else
        xbins=linspace(min_big2big_relative_distance,max_small2big_relative_distance,100);
    end
    [number_of_points_in_the_bin_small2big,~] = hist(small2big_relative_distance_vct,xbins);
    [number_of_points_in_the_bin_big2big,~] = hist(big2big_relative_distance_vct,xbins);
    discrete_probability_small2big=number_of_points_in_the_bin_small2big/sum(number_of_points_in_the_bin_small2big);
    discrete_probability_big2big=number_of_points_in_the_bin_big2big/sum(number_of_points_in_the_bin_big2big);
    Hellinger_distance(i)=(1/sqrt(2))*sqrt(sum((sqrt(discrete_probability_big2big)-sqrt(discrete_probability_small2big)).^2));
    if (abs(min_small2big_relative_distance)<10 && length(field{parallel_plane(i,2)}.points_of_filed)> (r*c/13))||(Hellinger_distance(i) < threshold_align_plane)
        aligned_plane(times_of_alignment,:)=[parallel_plane(i,1),parallel_plane(i,2)];
        aligned_plane_relative_distance(times_of_alignment,1)=small2big_sum_of_relative_distance;
        times_of_alignment=times_of_alignment+1;
    end
end

aligned_plane=unique(aligned_plane,'rows');  %delete the repeated rows
return%|||||||||||||||||||||||||||||||||||||||||||||||||||

function [D_square_vct,Relative_distance_all_points,Relative_distance_each_point ] = measure_points_distances_from_plane(points, equations_of_plan)

[A] = equations_of_plan(1);
[B] = equations_of_plan(2);
[C] = equations_of_plan(3);
[D] = equations_of_plan(4);

D_square_vct = (points * [A B C]' + D).^2 / ([A B C] *[A B C]'); %TT
Relative_distance_all_points = sum(points * [A B C]' + D) / sqrt([A B C] *[A B C]');
Relative_distance_each_point = (points * [A B C]' + D) / sqrt([A B C] *[A B C]');
return %||||||||||||||||||||||||||||||||||||||||||||||||||

function [field_index,merged_plane,times_of_merge,field_delete_merged_plane]=merge_surface(field,aligned_plane,field_index)
merged_plane=[];
times_of_merge=1;
se=strel('square',2);
field_delete_merged_plane=field;
for   i=1: length(aligned_plane)
    if ~(isempty(field{aligned_plane(i,1)}.points_of_filed) && isempty(field{aligned_plane(i,2)}.points_of_filed))
        ang_of_two_planes=calculate_angle(field{aligned_plane(i,1)}.equations_of_plan_{end}(1:4),field{aligned_plane(i,2)}.equations_of_plan_{end}(1:4));
        if ang_of_two_planes <25 %angle=25 for scene1 only
            mask1=field_index==aligned_plane(i,1);
            mask2=field_index==aligned_plane(i,2);
            mask1_dilate=imdilate(mask1,se);
            mask2_dilate=imdilate(mask2,se);
            if sum((mask1_dilate(:)&mask2_dilate(:)))~=0 % the planes share the same points?
                merged_plane(times_of_merge,:)=[aligned_plane(i,1),aligned_plane(i,2)];
                times_of_merge=times_of_merge+1;
                %maintain the plane which has more points
                new_plane_idx=aligned_plane(i,2);
                field_delete_merged_plane{1, new_plane_idx}.points_of_filed=vertcat(field{1, aligned_plane(i,1)}.points_of_filed, field{1, aligned_plane(i,2)}.points_of_filed);
                ri_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 1);
                ci_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 2);
                z_data = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 3);
                sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');
                A =  sf.p10;
                B = sf.p01;
                C = -1;
                D = sf.p00;
                field_delete_merged_plane{1, new_plane_idx}.equations_of_merged_plane = [A, B, C, D];
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                field_delete_merged_plane{1, aligned_plane(i,1)}.points_of_filed=[];
                field_delete_merged_plane{1, aligned_plane(i,1)}.equations_of_plan_{length(field{1, aligned_plane(i,1)}.equations_of_plan_)}=[nan,nan,nan,nan,nan];
                [r1,c1]=find(field_index==aligned_plane(i,1));
                field_index(sub2ind(size(field_index), r1, c1))=new_plane_idx; %the index of be merged plane is absent from this matrix. e.g is plane 2 has been merged into plane 1, then the original index 2 will be changed to 1. And no index value 2 exist any more.
                aligned_plane(i,:)=0;
                %aligned_plane(aligned_plane==aligned_plane(i,1))=aligned_plane(i,2);
            end
        end
    end
end
for   i=1: length(aligned_plane)
    if aligned_plane(i,1)~=0
        ang_of_two_planes=calculate_angle(field{aligned_plane(i,1)}.equations_of_plan_{end}(1:4),field{aligned_plane(i,2)}.equations_of_plan_{end}(1:4));
        if ang_of_two_planes <33
            mask1=field_index==aligned_plane(i,1);
            mask2=field_index==aligned_plane(i,2);
            mask1_dilate=imdilate(mask1,se);
            mask2_dilate=imdilate(mask2,se);
            if sum((mask1_dilate(:)&mask2_dilate(:)))~=0 % the planes share the same points?
                merged_plane(times_of_merge,:)=[aligned_plane(i,1),aligned_plane(i,2)];
                times_of_merge=times_of_merge+1;
                maintain the plane which has more points
                new_plane_idx=aligned_plane(i,2);
                field_delete_merged_plane{1, new_plane_idx}.points_of_filed=vertcat(field{1, aligned_plane(i,1)}.points_of_filed, field{1, aligned_plane(i,2)}.points_of_filed);
                ri_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 1);
                ci_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 2);
                z_data = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 3);
                sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');
                A =  sf.p10;
                B = sf.p01;
                C = -1;
                D = sf.p00;
                field_delete_merged_plane{1, new_plane_idx}.equations_of_merged_plane = [A, B, C, D];
                field_delete_merged_plane{1, new_plane_idx}.equations_of_merged_plane=(field{aligned_plane(i,1)}.Normal_vector_final(1:end-1) + field{aligned_plane(i,2)}.Normal_vector_final(1:end-1))./2;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                field_delete_merged_plane{1, aligned_plane(i,1)}.points_of_filed=[];
                field_delete_merged_plane{1, aligned_plane(i,1)}.equations_of_plan_{length(field{1, aligned_plane(i,1)}.equations_of_plan_)}=[nan,nan,nan,nan,nan];
                [r1,c1]=find(field_index==aligned_plane(i,1));
                field_index(sub2ind(size(field_index), r1, c1))=new_plane_idx; %the index of be merged plane is absent from this matrix. e.g is plane 2 has been merged into plane 1, then the original index 2 will be changed to 1. And no index value 2 exist any more.
                aligned_plane(i,:)=0;
                aligned_plane(aligned_plane==aligned_plane(i,1))=aligned_plane(i,2);
            end
        end
    end
end
% 
% for   i=1: length(aligned_plane)
%     if aligned_plane(i,1)~=0
%         ang_of_two_planes=calculate_angle(field{aligned_plane(i,1)}.equations_of_plan_{end}(1:4),field{aligned_plane(i,2)}.equations_of_plan_{end}(1:4));
%         if ang_of_two_planes <11
%             mask1=field_index==aligned_plane(i,1);
%             mask2=field_index==aligned_plane(i,2);
%             mask1_dilate=imdilate(mask1,se);
%             mask2_dilate=imdilate(mask2,se);
%             if sum((mask1_dilate(:)&mask2_dilate(:)))~=0 % the planes share the same points?
%                 merged_plane(times_of_merge,:)=[aligned_plane(i,1),aligned_plane(i,2)];
%                 times_of_merge=times_of_merge+1;
%                 %maintain the plane which has more points
%                 new_plane_idx=aligned_plane(i,2);
%                 field_delete_merged_plane{1, new_plane_idx}.points_of_filed=vertcat(field{1, aligned_plane(i,1)}.points_of_filed, field{1, aligned_plane(i,2)}.points_of_filed);
%                 %%%%%%%%%%%%%%%% JZ modified on 20170113 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 ri_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 1);
%                 ci_idx = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 2);
%                 z_data = field_delete_merged_plane{1, new_plane_idx}.points_of_filed(:, 3);
%                 sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');
%                 A =  sf.p10;
%                 B = sf.p01;
%                 C = -1;
%                 D = sf.p00;
%                 field_delete_merged_plane{1, new_plane_idx}.equations_of_merged_plane = [A, B, C, D];
%                 % field_delete_merged_plane{1, new_plane_idx}.equations_of_merged_plane=(field{aligned_plane(i,1)}.Normal_vector_final(1:end-1) + field{aligned_plane(i,2)}.Normal_vector_final(1:end-1))./2;
%                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 field_delete_merged_plane{1, aligned_plane(i,1)}.points_of_filed=[];
%                 field_delete_merged_plane{1, aligned_plane(i,1)}.equations_of_plan_{length(field{1, aligned_plane(i,1)}.equations_of_plan_)}=[nan,nan,nan,nan,nan];
%                 [r1,c1]=find(field_index==aligned_plane(i,1));
%                 field_index(sub2ind(size(field_index), r1, c1))=new_plane_idx; %the index of be merged plane is absent from this matrix. e.g is plane 2 has been merged into plane 1, then the original index 2 will be changed to 1. And no index value 2 exist any more.
%                 aligned_plane(i,:)=0;
%                 %aligned_plane(aligned_plane==aligned_plane(i,1))=aligned_plane(i,2);
%             end
%         end
%     end
% end
return%|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

