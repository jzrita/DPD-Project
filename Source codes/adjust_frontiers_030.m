function [fields_data] = adjust_frontiers_030(fields_data_in, basic_penetrating_element_shape, epoch_idx, show)

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

fields_data = fields_data_in;
fields_data.basic_data.maximum_allowed_dissimilarity=3;
fields_data.historical_field_in_epoch_{epoch_idx}.epoch_idx = epoch_idx;
fields_data.historical_field_in_epoch_{epoch_idx}.field_ = fields_data.field_;
fields_data.historical_field_in_epoch_{epoch_idx}.fields = fields_data.fields;
%fields_data.historical_field_in_epoch_{epoch_idx}.fields = fields_data.parallel_surface_detection.field_index;
%fields_data.fields=fields_data.parallel_surface_detection.field_index;
%shared_elements_ = double.empty(0,0);
fields_mutual_data = [];

fields_data.flooting_fields.current_field = zeros(size(fields_data.fields));
fields_data.flooting_fields.previous_field = zeros(size(fields_data.fields));
fields_data.flooting_fields.u_field_caused_change = zeros(size(fields_data.fields));

idx_new_fields = max(max(fields_data.fields)) +1;
num_original_fields = length(fields_data.field_);

for i_field = 1 : num_original_fields
    [fields_mutual_data] = update_fields_mutual_data(fields_data, i_field, fields_mutual_data, num_original_fields);
    fields_data.field_{i_field}.field_type = 'original_field'; % just to destinguish the original field from the newly created one
    for u_field = 1 : num_original_fields
        shared_elements = double.empty(0,0);
        fields_mutual_data{i_field, u_field}.shared_elements_ = shared_elements;
        disp(sprintf('TTILLO: data i-field=%d; u-field=%d', i_field, u_field));
        if (fields_mutual_data{i_field, u_field}.neighbours == 1)
            if (fields_mutual_data{i_field, u_field}.may_need_process == 1)
                basic_penetrating_element.shape = basic_penetrating_element_shape;
                if (i_field == 10) && (u_field == 14)
                    % for debugging
                    dummy_var = 1;
                end;
                if (i_field == 10) && (u_field == 21)
                    % for debugging
                    dummy_var = 1;
                end;
                [basic_penetrating_element.radius_of_penetrating_element] = get_radius_of_penetrating_element(fields_data, u_field, i_field);
                
                %            disp('TTILLO: the shift_from_boarder is used here to ensure that the intersection lines (where the penetrating_element moves) is completly within the image');
                shift_from_boarder = basic_penetrating_element.radius_of_penetrating_element;
                [intersection_line_tmp, mask_of_intersection] = intersection_of_two_surfaces_010(fields_data, u_field, i_field, shift_from_boarder);
                if ~isnan(mask_of_intersection)
                fields_mutual_data{i_field, u_field}.intersection_line = intersection_line_tmp;
                DEBUG = 1;
                if DEBUG
                    if fields_mutual_data{i_field, u_field}.intersection_line.valid_intersection
                        figure(111); image(mod(fields_data.fields.*(~mask_of_intersection) + 30*mask_of_intersection, size(settelments_colormap, 1))); colormap(settelments_colormap); title(sprintf('The intersection line between i-field=%d and u-field=%d; radius-of-penetrating-element=%d', i_field, u_field, basic_penetrating_element.radius_of_penetrating_element));
%                    else
%                        figure(111); image(mod(fields_data.fields.*(~mask_of_intersection) + 30*mask_of_intersection, size(settelments_colormap, 1))); title(sprintf('Strange behaviour of i-field=%d and u-field=%d; radius-of-penetrating-element=%d', i_field, u_field, basic_penetrating_element.radius_of_penetrating_element));
                    end;
                end
                
                if fields_mutual_data{i_field, u_field}.intersection_line.valid_intersection
                    [basic_penetrating_element] = get_basic_penetrating_element(basic_penetrating_element, fields_mutual_data{i_field, u_field}.intersection_line, mask_of_intersection);
                    
                    [mask_of_penetrating_element_operation] = get_mask_of_penetrating_element_operation(fields_data, u_field, i_field, basic_penetrating_element, shift_from_boarder, show);
                    
                    order_of_checking_points = [1: size(fields_mutual_data{i_field, u_field}.intersection_line.points, 1)];
                    
                    [shared_elements_] = get_shared_structures(fields_data, u_field, i_field, fields_mutual_data{i_field, u_field}.intersection_line, basic_penetrating_element, show, order_of_checking_points);
                    fields_mutual_data{i_field, u_field}.shared_elements_ = shared_elements_;
                    %%%%% not here this defination
                    epoch_idx = 1;
                    attack_idx = 1;
                    
                    fields_mutual_data{i_field, u_field}.basic_penetrating_element = basic_penetrating_element;
                    
                    [fields_data, fields_mutual_data] = check_where_shared_elements_should_be(fields_data, fields_mutual_data, i_field, u_field, epoch_idx, show, attack_idx);
                    [fields_data, fields_mutual_data, move_points_from_side_to_side_has_occured] = move_points_from_side_to_side(fields_data, fields_mutual_data, i_field, u_field, epoch_idx, show, attack_idx);

                    if move_points_from_side_to_side_has_occured
                       [fields_data, fields_mutual_data, idx_new_fields] = check_if_any_isolated_groups(fields_data, i_field, u_field, fields_mutual_data, epoch_idx, show, attack_idx, idx_new_fields, mask_of_penetrating_element_operation);
                    end;
                else
%                    warning(sprintf('TTILLO: Strange behaviour of i-field=%d and u-field=%d : NO valid intersection line', i_field, u_field));
                    disp('This phynomena to be discribed in the paper')
                end;
                clear basic_penetrating_element;
            else
                warning(sprintf('TTILLO: no need for process of i-field=%d; u-field=%d, becasue already processed', i_field, u_field));    
                end
            end;
        end;
    end;
end;

fields_data.fields_mutual_data = fields_mutual_data;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [fields_mutual_data] = update_fields_mutual_data(fields_data, i_field, fields_mutual_data, num_original_fields)

if i_field == 1,
    fields_mutual_data = cell(num_original_fields, num_original_fields);
    for ii = 1 : num_original_fields
        for uu = 1 : num_original_fields
            fields_mutual_data{ii, uu}.neighbours = 0;
            fields_mutual_data{ii, uu}.may_need_process = 1;
        end;
    end;
end;

[is_island, neighbours_list, neighbours_map] = fields_touching_i_filed(fields_data, i_field);
    
length(neighbours_list)
for idx = 1 : length(neighbours_list)
    fields_mutual_data{i_field, neighbours_list(idx)}.neighbours = 1;
end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



function [fields_data, fields_mutual_data, idx_new_fields] = check_if_any_isolated_groups(fields_data, i_field, u_field, fields_mutual_data, epoch_idx, show, attack_idx, idx_new_fields, mask_of_penetrating_element_operation);

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

[r,c] = size(fields_data.fields);

if (i_field == 37) && (u_field == 7)    
    % for debugging
    dummy_var = 1;
end;

mask_of_initial_i_field = (fields_data.fields == i_field);

DEBUG = 0   ;
if DEBUG 
    figure; image(mask_of_initial_i_field*255); title('mask_of_initial_i_field');
end;

CC = bwconncomp(mask_of_initial_i_field);
mask_of_initial_i_field_excluding_main_i_field = mask_of_initial_i_field;
            
if CC.NumObjects > 1
    for jj = 1 : CC.NumObjects
        vct_tmp(jj) = length(CC.PixelIdxList{jj});
    end;
    [D,I] = max(vct_tmp);
    clear vct_tmp;
    mask_of_initial_i_field_excluding_main_i_field(CC.PixelIdxList{I}) = 0;
    if DEBUG 
        figure; image(mask_of_initial_i_field_excluding_main_i_field*255); title('mask_of_initial_i_field_excluding_main_i_field');
    end;
    CC1 = bwconncomp(mask_of_initial_i_field_excluding_main_i_field);
    for ii = 1 : CC1.NumObjects
        [fields_data, idx_new_fields] = re_allocate_isolated_groups(fields_data, CC1.PixelIdxList{ii}, i_field, u_field, epoch_idx, show, attack_idx, idx_new_fields);
    end;
end;

% this here to check fro single pixels width objects and try to eliminate
% them
% 
% mask_of_isolated_1 = bwhitmiss(mask_of_initial_i_field, [0 1 0], [1 0 1]);
% mask_of_isolated_2 = bwhitmiss(mask_of_initial_i_field, [0 1 0]', [1 0 1]');
% mask_of_isolated_3 = bwhitmiss(mask_of_initial_i_field, [0 1 ; 1 0], [1 0 ; 0 1]);
% mask_of_isolated_4 = bwhitmiss(mask_of_initial_i_field, [1 0 ; 0 1], [0 1 ; 1 0]);
% 
% mask_of_isolated = mask_of_isolated_1 | mask_of_isolated_2 | mask_of_isolated_3 | mask_of_isolated_4; 
% 
% mask_of_isolated = zeros(r,c);
% 
% CC = bwconncomp(mask_of_initial_i_field.* (~mask_of_isolated));
% mask_of_initial_i_field_excluding_main_i_field = mask_of_initial_i_field;
% 
% if (i_field == 37) && (u_field == 7)    
%     % for debugging
%     dummy_var = 1;
% end;
%             
% if CC.NumObjects > 1
%     for jj = 1 : CC.NumObjects
%         vct_tmp(jj) = length(CC.PixelIdxList{jj});
%     end;
%     [D,I] = max(vct_tmp);
%     clear vct_tmp;
%     mask_of_initial_i_field_excluding_main_i_field(CC.PixelIdxList{I}) = 0;
%     if DEBUG 
%         figure; image(mask_of_initial_i_field_excluding_main_i_field*255); title('mask_of_initial_i_field_excluding_main_i_field');
%     end;
%     CC1 = bwconncomp(mask_of_initial_i_field_excluding_main_i_field.* (~mask_of_isolated).* mask_of_penetrating_element_operation);
%     if DEBUG 
%         figure; image(mask_of_initial_i_field_excluding_main_i_field.* (~mask_of_isolated).* mask_of_penetrating_element_operation*255); title('mask_of_initial_i_field_excluding_main_i_field.* (~mask_of_isolated).* mask_of_penetrating_element_operation');
%     end;
%     for ii = 1 : CC1.NumObjects
%         [fields_data, idx_new_fields] = re_allocate_isolated_groups(fields_data, CC1.PixelIdxList{ii}, i_field, u_field, epoch_idx, show, attack_idx, idx_new_fields);
%     end;
% end;
% 
% CC2 = bwconncomp(mask_of_isolated.* mask_of_penetrating_element_operation);
% for ii = 1 : CC2.NumObjects
%     [fields_data, idx_new_fields] = re_allocate_isolated_groups(fields_data, CC2.PixelIdxList{ii}, i_field, u_field, epoch_idx, show, attack_idx, idx_new_fields);
% end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [mask_of_operation] = get_mask_of_penetrating_element_operation(fields_data, u_field, i_field, basic_penetrating_element, shift_from_boarder, show)

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

[r,c] = size(fields_data.fields);

[intersection_line_tmp, full_mask_of_intersection] = intersection_of_two_surfaces_010(fields_data, u_field, i_field, 0);

element{1} = basic_penetrating_element.element1;
element{2} = basic_penetrating_element.element2;
element{3} = basic_penetrating_element.element3;

mask_of_operation = zeros(r + 4*shift_from_boarder, c + 4*shift_from_boarder);
mask_of_operation(2*shift_from_boarder + [1:r], 2*shift_from_boarder + [1:c]) = full_mask_of_intersection;

final_mask = zeros(r + 4*shift_from_boarder, c + 4*shift_from_boarder);
for ii = 1 : 3
    for pii = 1 : size(element{ii}, 1)
        mask_tmp = zeros(r + 4*shift_from_boarder, c + 4*shift_from_boarder);
        
        mask_tmp(2*shift_from_boarder + [1:r] + element{ii}(pii, 1), 2*shift_from_boarder + [1:c] + element{ii}(pii, 2)) = full_mask_of_intersection;
        final_mask = final_mask + mask_tmp;
    end;
end;

mask_of_operation = final_mask(2*shift_from_boarder + [1:r], 2*shift_from_boarder + [1:c]);
mask_of_operation = bwmorph(mask_of_operation, 'fill');

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [basic_penetrating_element] = get_basic_penetrating_element(basic_penetrating_element, intersection_line, mask_of_intersection);
% Be careful x and y here are swapped :-(
radius_of_penetrating_element = basic_penetrating_element.radius_of_penetrating_element;

switch basic_penetrating_element.shape
    case {'disk'}
        xx = [1: 2*radius_of_penetrating_element+1];
        xx = repmat(xx, 2*radius_of_penetrating_element+1, 1);
        [r,c] = size(xx);
        
        xx = xx - radius_of_penetrating_element -1;
        
        xx(1,1) = inf;
        xx(1,c) = inf;
        xx(r,1) = inf;
        xx(r,c) = inf;
        
        yy = xx';
        
        xx = xx(:);
        yy = yy(:);
        
        idx_tmp = xx ~= inf;
        xx = xx(idx_tmp);
        yy = yy(idx_tmp);
        
        idx_tmp = yy ~= inf;
        xx = xx(idx_tmp);
        yy = yy(idx_tmp);
        
        basic_penetrating_element.whole(:, 2) = xx;
        basic_penetrating_element.whole(:, 1) = yy;
        
        [r,c] = size(mask_of_intersection);
        %mask = zeros(r,c);
        
        theta = atan2d(yy, xx);
        theta = 360*(theta < 0) + theta;
        
        rr_lin1 = atan2d(intersection_line.tri_d_line_equation.rr(1), intersection_line.tri_d_line_equation.rr(2));
        rr_lin1 = 360*(rr_lin1 < 0) + rr_lin1;
        
        rr_lin2 = atan2d(-intersection_line.tri_d_line_equation.rr(1), -intersection_line.tri_d_line_equation.rr(2));
        rr_lin2 = 360*(rr_lin2 < 0) + rr_lin2;
        
        idx_tmp = (min(rr_lin1, rr_lin2) < theta) & (theta < max(rr_lin1, rr_lin2));
        
        basic_penetrating_element.element1(:, 2) = xx(~idx_tmp);
        basic_penetrating_element.element1(:, 1) = yy(~idx_tmp);
        %basic_penetrating_element.element1.cc = xx(~idx_tmp);
        %basic_penetrating_element.element1.rr = yy(~idx_tmp);
        
        
        idx_tmp = (max(rr_lin1, rr_lin2) < theta) | (theta < min(rr_lin1, rr_lin2));
        
        basic_penetrating_element.element2(:, 2) = xx(~idx_tmp);
        basic_penetrating_element.element2(:, 1) = yy(~idx_tmp);
        %basic_penetrating_element.element2.cc = xx(~idx_tmp);
        %basic_penetrating_element.element2.rr = yy(~idx_tmp);
        
        if abs(intersection_line.tri_d_line_equation.rr(1)) <= abs(intersection_line.tri_d_line_equation.rr(2))
            slop = intersection_line.tri_d_line_equation.rr(1)/ intersection_line.tri_d_line_equation.rr(2);
            intersection_line(:, 2) = [-radius_of_penetrating_element : radius_of_penetrating_element];
            intersection_line(:, 1) = round(slop * (intersection_line(:, 2)) );
        else
            % to be checked
            error('TTillo: Please debug this part')
            slop = intersection_line.tri_d_line_equation.rr(2)/ intersection_line.tri_d_line_equation.rr(1);
            intersection_line(:, 1) = [-radius_of_penetrating_element : radius_of_penetrating_element];
            intersection_line(:, 2) = round(slop * (intersection_line(:, 1)) );
        end;
        
        basic_penetrating_element.element1 = union(basic_penetrating_element.element1, intersection_line, 'rows');
        basic_penetrating_element.element2 = union(basic_penetrating_element.element2, intersection_line, 'rows');
        
    case {'line'}
        %         xx = [1: 2*radius_of_penetrating_element+1];
        %         xx = repmat(xx, 2*radius_of_penetrating_element+1, 1);
        %         [r,c] = size(xx);
        %
        %         xx = xx - radius_of_penetrating_element -1;
        %         yy = xx';
        %
        %         % perpendacular on the penetrating line
        
        slope_of_line = -intersection_line.tri_d_line_equation.rr(2)/ intersection_line.tri_d_line_equation.rr(1);
        
        if abs(slope_of_line) <= 1
            %  error('TTillo: Please debug this part')
            basic_penetrating_element.element1(:, 2) = [-radius_of_penetrating_element : -1];
            basic_penetrating_element.element1(:, 1) = round(slope_of_line * (basic_penetrating_element.element1(:, 2)) );
            basic_penetrating_element.element2(:, 2) = [1 : radius_of_penetrating_element];
            basic_penetrating_element.element2(:, 1) = round(slope_of_line * (basic_penetrating_element.element2(:, 2)) );
            basic_penetrating_element.element3(:, 2) = [0];
            basic_penetrating_element.element3(:, 1) = round(slope_of_line * (basic_penetrating_element.element3(:, 2)) );
        else
            % to be checked
            %error('TTillo: Please debug this part')
            basic_penetrating_element.element1(:, 1) = [-radius_of_penetrating_element : -1];
            basic_penetrating_element.element1(:, 2) = round((1/slope_of_line) * (basic_penetrating_element.element1(:, 1)) );
            basic_penetrating_element.element2(:, 1) = [1 : radius_of_penetrating_element];
            basic_penetrating_element.element2(:, 2) = round((1/slope_of_line) * (basic_penetrating_element.element2(:, 1)) );
            basic_penetrating_element.element3(:, 1) = [0];
            basic_penetrating_element.element3(:, 2) = round((1/slope_of_line) * (basic_penetrating_element.element3(:, 1)) );
        end;
        
        basic_penetrating_element.whole = union(basic_penetrating_element.element1, basic_penetrating_element.element2, 'rows');
        basic_penetrating_element.whole = union(basic_penetrating_element.whole, basic_penetrating_element.element3, 'rows');
        
    otherwise
        error('TTILLO: non supported shape of penetrating element')
end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [radius_of_penetrating_element] = get_radius_of_penetrating_element(fields_data, u_field, i_field);

equations_of_i_field = get_equations_of_plan(fields_data, i_field);
equations_of_u_field = get_equations_of_plan(fields_data, u_field);

[A] = equations_of_i_field(1);
[B] = equations_of_i_field(2);
[C] = equations_of_i_field(3);
[D] = equations_of_i_field(4);
n1 = [A B C];
n1 = n1 / norm(n1, 2);
d1 = -D;

[A] = equations_of_u_field(1);
[B] = equations_of_u_field(2);
[C] = equations_of_u_field(3);
[D] = equations_of_u_field(4);
n2 = [A B C];
n2 = n2 / norm(n2, 2);
d2 = -D;

theta_between_plans = acosd( dot(n1, n2) );

%radius_of_penetrating_element = 6;
radius_of_penetrating_element = ceil(abs(fields_data.basic_data.maximum_allowed_dissimilarity / tand(theta_between_plans)));
radius_of_penetrating_element = min(radius_of_penetrating_element +3, 12);
%radius_of_penetrating_element = min(radius_of_penetrating_element +6, 20);
warning('TTILLO, why 12, in the above euation')

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [shared_elements_] = get_shared_structures(fields_data, u_field, i_field, intersection_line, basic_penetrating_element, show, order_of_checking_points)

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

[r,c] = size(fields_data.fields);

basic_penetrating_element_tmp{1} = basic_penetrating_element.element1;
basic_penetrating_element_tmp{2} = basic_penetrating_element.element2;
basic_penetrating_element_tmp{3} = basic_penetrating_element.element3; % central point
             
mask_of_elements = zeros(r, c);
mask_of_elements_type = zeros(r, c);
element_idx = 0;

for ii = 1 : 2
    
    already_in_fields = 'NO';

    for jj = 1 : length(order_of_checking_points)
        
        pp = order_of_checking_points(jj);
        
        ponits_tmp(:, 1) = (basic_penetrating_element_tmp{ii}(:,1) + intersection_line.points(pp, 1));
        ponits_tmp(:, 2) = (basic_penetrating_element_tmp{ii}(:,2) + intersection_line.points(pp, 2));

        idx_tmp = sub2ind([r, c], ponits_tmp(:, 1), ponits_tmp(:, 2));
        ponits_tmp(:, 4) = fields_data.fields(idx_tmp);
        ponits_tmp(:, 3) = fields_data.im_d(idx_tmp);
        clear idx_tmp;

        DEBUG = 0;
        if DEBUG
            mask = zeros(r, c);
            mask(sub2ind([r, c], ponits_tmp(:, 1), ponits_tmp(:, 2))) = 1;
            figure(90); image(255*mask); colormap(gray(2)); 
            clear mask;
        end;
        
        if any(ponits_tmp(:, 4) == i_field) && any(ponits_tmp(:, 4) == u_field)
            in_fields = 'IU';
        elseif any(ponits_tmp(:, 4) == i_field) && any(ponits_tmp(:, 4) ~= u_field)
            in_fields = 'Ix';
        elseif any(ponits_tmp(:, 4) ~= i_field) && any(ponits_tmp(:, 4) == u_field)
            in_fields = 'xU';
        else any(ponits_tmp(:, 4) ~= i_field) && any(ponits_tmp(:, 4) ~= u_field);
            in_fields = 'xx';
        end;       
               
        if all(already_in_fields(end+[-1:0]) == in_fields)
            shared_elements_{element_idx}.element_length = shared_elements_{element_idx}.element_length + 1;
 
            ponits_tmp = ponits_tmp(ponits_tmp(:, 4) == i_field, :);
            shared_elements_{element_idx}.points = union(shared_elements_{element_idx}.points, ponits_tmp, 'rows');

            switch basic_penetrating_element.shape
                case {'disk'}
                    error('TTILLO this case is not handled now; for the boarder of the penetrating objects');
                case {'line'}
                    if isempty(ponits_tmp)
                        ponits_tmp_of_boarder = ponits_tmp;
                    else
                        if ii == 1,
                            ponits_tmp_of_boarder = ponits_tmp(1, :);
                        elseif ii == 2,
                            ponits_tmp_of_boarder = ponits_tmp(end, :);
                        else
                            error('TTILLO this case is not handled now');
                        end;
                    end;
                    shared_elements_{element_idx}.points_of_boarder = union(shared_elements_{element_idx}.points_of_boarder, ponits_tmp_of_boarder, 'rows');
                    clear ponits_tmp_of_boarder;
            end;
            
            points_of_intersection_line_tmp = [intersection_line.points(pp, 1:3)];
            points_of_intersection_line_tmp(4) = fields_data.fields(points_of_intersection_line_tmp(1), points_of_intersection_line_tmp(2));
            shared_elements_{element_idx}.points_of_intersection_line = union(shared_elements_{element_idx}.points_of_intersection_line, points_of_intersection_line_tmp, 'rows');
        else
            element_idx = element_idx +1;

            already_in_fields = in_fields;
            shared_elements_{element_idx}.element_type = in_fields;

            shared_elements_{element_idx}.start = pp;            
            shared_elements_{element_idx}.element_length = 1;
 
            ponits_tmp = ponits_tmp(ponits_tmp(:, 4) == i_field, :);
            shared_elements_{element_idx}.points = ponits_tmp;
                        
            switch basic_penetrating_element.shape
                case {'disk'}
                    error('TTILLO this case is not handled now; for the boarder of the penetrating objects');
                case {'line'}
                    if isempty(ponits_tmp)
                        ponits_tmp_of_boarder = ponits_tmp;
                    else
                        if ii == 1,
                            ponits_tmp_of_boarder = ponits_tmp(1, :);
                        elseif ii == 2,
                            ponits_tmp_of_boarder = ponits_tmp(end, :);
                        else
                            error('TTILLO this case is not handled now');
                        end;
                    end;
                    shared_elements_{element_idx}.points_of_boarder = ponits_tmp_of_boarder;
                    clear ponits_tmp_of_boarder;
            end;
            
            points_of_intersection_line_tmp = [intersection_line.points(pp, 1:3)];
            points_of_intersection_line_tmp(4) = fields_data.fields(points_of_intersection_line_tmp(1), points_of_intersection_line_tmp(2));
            shared_elements_{element_idx}.points_of_intersection_line = points_of_intersection_line_tmp;
            
            shared_elements_{element_idx}.i_field = i_field;
            shared_elements_{element_idx}.u_field = u_field;
            
        end;
        DEBUG = 0;
        if DEBUG
            mask_of_elements(sub2ind([r, c], shared_elements_{element_idx}.points(:, 1), shared_elements_{element_idx}.points(:, 2))) = element_idx;
            figure(91); image(mask_of_elements); colormap(settelments_colormap); title(sprintf('The current shared elements=%d', element_idx));
            switch shared_elements_{element_idx}.element_type
                case 'IU',
                    mask_of_elements_type(sub2ind([r, c], shared_elements_{element_idx}.points(:, 1), shared_elements_{element_idx}.points(:, 2))) = 16;                    
                case 'Ix'
                    mask_of_elements_type(sub2ind([r, c], shared_elements_{element_idx}.points(:, 1), shared_elements_{element_idx}.points(:, 2))) = 8;
                case 'xU',
                    mask_of_elements_type(sub2ind([r, c], shared_elements_{element_idx}.points(:, 1), shared_elements_{element_idx}.points(:, 2))) = 4;
                case 'xx',
                    mask_of_elements_type(sub2ind([r, c], shared_elements_{element_idx}.points(:, 1), shared_elements_{element_idx}.points(:, 2))) = 2;
            end;
            figure(92); image(mask_of_elements_type); colormap(settelments_colormap); title(sprintf('The mask_of_elements_type'));            
        end;
        clear ponits_tmp
        clear points_of_intersection_line_tmp
    end;
    
end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [is_island, neighbours_list, neighbours_map] = fields_touching_i_filed(fields_data, i_field)

mask_of_current_batch = (fields_data.fields == i_field);
boarder_of_current_field = xor(mask_of_current_batch, imdilate(mask_of_current_batch, ones(3)));
neighbours_map = fields_data.fields.* boarder_of_current_field;

neighbours_list = unique(neighbours_map(:));

neighbours_list = neighbours_list(neighbours_list ~= 0);

is_island = isempty( neighbours_list );

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [fields_data, fields_mutual_data, move_points_from_side_to_side_has_occured] = move_points_from_side_to_side(fields_data, fields_mutual_data, i_field, u_field, epoch_idx, show, attack_idx)

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

shared_elements_tmp = fields_mutual_data{i_field, u_field}.shared_elements_;

[r,c] = size(fields_data.fields);

move_points_from_side_to_side_has_occured = 0;

for ii = 1 : length(shared_elements_tmp)            
    
    if fields_mutual_data{i_field, u_field}.shared_elements_{ii}.move_from_I_to_U
        
        fields_mutual_data{u_field, i_field}.may_need_process = 0;        
        move_points_from_side_to_side_has_occured = 1;

        points = fields_mutual_data{i_field, u_field}.shared_elements_{ii}.points(:, 1:3);        
        points = vertcat(fields_mutual_data{i_field, u_field}.shared_elements_{ii}.points_of_intersection_line(:, 1:3), points);
%        t1 = size(points, 1);
%        idx = (fields_data.fields == u_field);
%        t4 = sum(sum(idx));

        idx_tmp = sub2ind([r,c], points(:, 1), points(:, 2));
        fields_data.fields(idx_tmp) = u_field;
%        t2 = size(fields_data.field_{u_field}.points_of_filed, 1);

        idx = (fields_data.fields == u_field);
        [points_of_filed_tmp(:, 1), points_of_filed_tmp(:, 2)] = find(idx);
        points_of_filed_tmp(:, 3) = fields_data.im_d( idx );    
        fields_data.field_{u_field}.points_of_filed = points_of_filed_tmp;
%       t3 = size(fields_data.field_{u_field}.points_of_filed, 1);
        clear points_of_filed_tmp
        
%        disp(sprintf('t1=%d; t2=%d (t4=%d); t1+t2 = %d; %d', t1, t2, t4, t1+ t2, t3 ))
        
        idx = (fields_data.fields == i_field);
        [points_of_filed_tmp(:, 1), points_of_filed_tmp(:, 2)] = find(idx);
        points_of_filed_tmp(:, 3) = fields_data.im_d( idx );
        fields_data.field_{i_field}.points_of_filed = points_of_filed_tmp;
        clear points_of_filed_tmp

%         clear i_field_points;
%         if 0, % refine_plan_parameters;
%             warning('TTILLO: here we may update the plan''s equation');
%         end;

points_of_filed_for_update = fields_data.field_{i_field}.points_of_filed;

equations_of_plan_for_update = fields_data.field_{i_field}.equations_of_plan_;
bridth_of_grow = equations_of_plan_for_update{ length(equations_of_plan_for_update)}(5);

[fields_data.field_{i_field}.equations_of_plan_] = refine_plan_parameters(equations_of_plan_for_update, points_of_filed_for_update, bridth_of_grow);
   
        DEBUG = 0;
        if DEBUG
            figure(40); image(fields_data.fields); colormap(settelments_colormap); title(sprintf('epoch=%d : move from %d to %d', epoch_idx, u_field, i_field));
        end
        
        attack_idx = attack_idx +1;
    end;
end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [equations_of_plan_] = refine_plan_parameters(equations_of_plan_, points_of_filed, bridth_of_grow)

[A] = equations_of_plan_{ length(equations_of_plan_) }(1);
[B] = equations_of_plan_{ length(equations_of_plan_) }(2);
[C] = equations_of_plan_{ length(equations_of_plan_) }(3);
[D] = equations_of_plan_{ length(equations_of_plan_) }(4);
if size(points_of_filed,1)>=3
%    d = -ax_0-by_0-cz_0
ri_idx = points_of_filed(:, 1);
ci_idx = points_of_filed(:, 2);
z_data = points_of_filed(:, 3);
sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');

% [ A, B, C, D] = [ sf.p10, sf.p01, -1, sf.p00];
%     [ sf.p10, sf.p01, -1, sf.p00]./[ A, B, C, D]
%     [ A, B, C, D]
d_new = sum((z_data - feval(sf, ri_idx, ci_idx)).^2);

sf_or = sf;
sf_or.p10 = A/(-C);
sf_or.p01 = B/(-C);
sf_or.p00 = D/(-C);
d_old = sum((z_data - feval(sf_or, ri_idx, ci_idx)).^2);

if d_new < d_old
    %        fprintf('bridth_of_grow is %d, num points %d, d_new=%.2f, d_old=%0.2f\n', bridth_of_grow, length(ri_idx), d_new, d_old);
    
    %      Linear model Poly11:
    %      sf(x,y) = p00 + p10*x + p01*y
    %      z = p00 + p10*x + p01*y
    %      p00 = - p10*x - p01*y + z
    A =  sf.p10;
    B = sf.p01;
    C = -1;
    D = sf.p00;
    equations_of_plan_{ length(equations_of_plan_) +1 } = [A, B, C, D, bridth_of_grow];
    %    [ sf.p10, sf.p01, -1, sf.p00]
end;
end
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [fields_data, idx_new_fields] = re_allocate_isolated_groups(fields_data, isolated_points_PixelIdxList, i_field, u_field, epoch_idx, show, attack_idx, idx_new_fields)
    
[r,c] = size(fields_data.im_d);
move_points_from_i_to_u_field_tmp = 0;

[isolated_points(:,1), isolated_points(:,2)] = ind2sub([r,c], isolated_points_PixelIdxList);
isolated_points(:,3) = fields_data.im_d(isolated_points_PixelIdxList);

% check if isolated_points are all englobed by u_field
if size(isolated_points, 1) <= 1    
    mask_of_current_batch = zeros(r,c);
    mask_of_current_batch( isolated_points_PixelIdxList ) = 1;    
    
    boarder_of_current_settelment = xor(mask_of_current_batch, imdilate(mask_of_current_batch, ones(3)));
    enemies_map = fields_data.fields.* boarder_of_current_settelment;
    
    enemies_list = unique(enemies_map(:));
    enemies_list = setdiff(enemies_list, 0);
    enemies_list = setdiff(enemies_list, i_field);
    
    if all(enemies_list == u_field)
        move_points_from_i_to_u_field_tmp = 1;
    end;
else
    equations_of_u_field = get_equations_of_plan(fields_data, u_field);
    [isolated_points_distance_from_u_field] = get_points_distance_from_plan(isolated_points(:,1:3), equations_of_u_field);
    
    [u_points_distance_from_u_field] = get_points_distance_from_plan(fields_data.field_{u_field}.points_of_filed, equations_of_u_field);
    
    move_points_from_i_to_u_field_tmp = mean(isolated_points_distance_from_u_field) <= max(u_points_distance_from_u_field);    
end;

if move_points_from_i_to_u_field_tmp,  % move from i into u_field
    fields_data.fields(isolated_points_PixelIdxList) = u_field;
    idx = (fields_data.fields == u_field);
    [points_of_filed_tmp(:, 1), points_of_filed_tmp(:, 2)] = find(idx);
    points_of_filed_tmp(:, 3) = fields_data.im_d( idx );    
    fields_data.field_{u_field}.points_of_filed = points_of_filed_tmp;

    fields_data.flooting_fields.current_field(isolated_points_PixelIdxList) = u_field;
    fields_data.flooting_fields.previous_field(isolated_points_PixelIdxList) = i_field;
    fields_data.flooting_fields.u_field_caused_change(isolated_points_PixelIdxList) = u_field;
    clear points_of_filed_tmp;
    
else % don't move i-(points) into u_field, but creat new field
    if size(isolated_points, 1) > 5
        warning('TTILLO: the minimum number of pixels of a newly created field is set here to 5')
        fields_data.fields(isolated_points_PixelIdxList) = idx_new_fields;
        fields_data.field_{idx_new_fields}.points_of_filed = isolated_points(:, 1:3);
        fields_data.field_{idx_new_fields}.field_type = 'flooting_field'; % i.e., created_after_adjust_frontier
        
        fields_data.flooting_fields.current_field(isolated_points_PixelIdxList) = idx_new_fields;
        fields_data.flooting_fields.previous_field(isolated_points_PixelIdxList) = i_field;
        fields_data.flooting_fields.u_field_caused_change(isolated_points_PixelIdxList) = u_field;
        [equations_of_plan_] = get_plan_equation(fields_data.field_{idx_new_fields}.points_of_filed);
        fields_data.field_{idx_new_fields}.equations_of_plan_{1} = equations_of_plan_;
        
    else % don't move i-(points) into u_field, but creat new zero-field (nothing)
        fields_data.fields(isolated_points_PixelIdxList) = 0;    
        fields_data.flooting_fields.current_field(isolated_points_PixelIdxList) = 0;
        fields_data.flooting_fields.previous_field(isolated_points_PixelIdxList) = i_field;
        fields_data.flooting_fields.u_field_caused_change(isolated_points_PixelIdxList) = u_field;
        
    end;
end;

% This here is to update i_field

idx = (fields_data.fields == i_field);
[points_of_filed_tmp(:, 1), points_of_filed_tmp(:, 2)] = find(idx);
points_of_filed_tmp(:, 3) = fields_data.im_d( idx );
fields_data.field_{i_field}.points_of_filed = points_of_filed_tmp;

idx_new_fields = idx_new_fields +1;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [equations_of_plan_] = get_plan_equation(points_of_filed)

%    d = -ax_0-by_0-cz_0
ri_idx = points_of_filed(:, 1);
ci_idx = points_of_filed(:, 2);
z_data = points_of_filed(:, 3);
sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');

A =  sf.p10;
B = sf.p01;
C = -1;
D = sf.p00;

bridth_of_grow = 1;

equations_of_plan_ = [A, B, C, D, bridth_of_grow];

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [fields_data, fields_mutual_data] = check_where_shared_elements_should_be(fields_data, fields_mutual_data, i_field, u_field, epoch_idx, show, attack_idx)
    
show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

shared_elements_tmp = fields_mutual_data{i_field, u_field}.shared_elements_;

for ii = 1 : length(shared_elements_tmp)            
    
    fields_mutual_data{i_field, u_field}.shared_elements_{ii}.move_from_I_to_U = 0;

    if all(shared_elements_tmp{ii}.element_type == 'IU')
        if size(shared_elements_tmp{ii}.points_of_boarder, 1) > 1
            
            equations_of_u_field = get_equations_of_plan(fields_data, u_field);
            [points_of_boarder_distance_from_u_field] = get_points_distance_from_plan(shared_elements_tmp{ii}.points_of_boarder(:,1:3), equations_of_u_field);
            
            equations_of_i_field = get_equations_of_plan(fields_data, i_field);
            [points_of_boarder_distance_from_i_field] = get_points_distance_from_plan(shared_elements_tmp{ii}.points_of_boarder(:,1:3), equations_of_i_field);

            [points_of_intersection_distance_from_u_field] = get_points_distance_from_plan(shared_elements_tmp{ii}.points_of_intersection_line(:,1:3), equations_of_u_field);
            [points_of_intersection_distance_from_i_field] = get_points_distance_from_plan(shared_elements_tmp{ii}.points_of_intersection_line(:,1:3), equations_of_i_field);
            DEBUG = 0;
            if DEBUG
                figure(30); 
                plot(points_of_boarder_distance_from_u_field,'r');
                hold on;
                plot(points_of_boarder_distance_from_i_field,'g');
                plot(points_of_intersection_distance_from_u_field,'.r'); 
                plot(points_of_intersection_distance_from_i_field,'.g'); 
                legend(strvcat('u-field boarder', 'i-field boarder', 'u-field center', 'i-field center')); 
                hold off;
            end;
            
            if mean(points_of_boarder_distance_from_i_field) > mean(points_of_boarder_distance_from_u_field)
                abs_diff_distance_of_boarder = mean(points_of_boarder_distance_from_i_field) - mean(points_of_boarder_distance_from_u_field);
                abs_diff_distance_of_center = abs( mean(points_of_intersection_distance_from_i_field) - mean(points_of_intersection_distance_from_u_field) );
                
                if abs_diff_distance_of_center < abs_diff_distance_of_boarder,                    
                    fields_mutual_data{i_field, u_field}.shared_elements_{ii}.move_from_I_to_U = 1;
                end;                
            end;            
        end;        
    end;
end;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [points_distance_from_plan] = get_points_distance_from_plan(points, equations_of_plan_);

points = points(:, 1:3);

[A] = equations_of_plan_(1);
[B] = equations_of_plan_(2);
[C] = equations_of_plan_(3);
[D] = equations_of_plan_(4);

D_square_vct = (points * [A B C]' + D).^2 / ([A B C] *[A B C]');

points_distance_from_plan = D_square_vct;

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



function [equations_of_plan_] = get_equations_of_plan(fields_data, current_settelment);

% for now we just use the most updated equation to represent the plan

equations_of_plan_ = fields_data.field_{current_settelment}.equations_of_plan_;
equations_of_plan_ = equations_of_plan_{length(equations_of_plan_)};

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

HDD = 'C:\Users\Tillo\Desktop\to del\1. Work to del' % Fujitsu
% HDD = 'I:\TT\1. Work'
%HDD = 'H:\TT\1. Work'

epoch_idx = 1;

%load(sprintf('%s/TT/1. Work/4. Matlab/DDPSD/results/from JinZhi/Indoor_S1-1_diamond_009_ascendorder_dissimilarity_2.mat', HDD), 'fields_data');
load(sprintf('%s/4. Matlab/DDPSD/results/from JinZhi/Indoor_S1-1_diamond_005_ascendorder.mat', HDD), 'fields_data');
load(sprintf('%s/4. Matlab/DDPSD/results/from JinZhi/stair2_square_016_ascendorder.mat', HDD), 'fields_data');
load(sprintf('%s/4. Matlab/DDPSD/results/from JinZhi/Indoor_S1-7_square_009_ascendorder_dissimilarity_3.mat', HDD), 'fields_data');


fields_data_in = fields_data;
basic_penetrating_element_shape = 'line';

show.show_im = 0;
show.fig_idx = 10;
load(sprintf('%s/4. Matlab/DDPSD/settelments_colormap.mat', HDD));
show.settelments_colormap = settelments_colormap;
figure(1); image(fields_data.fields); colormap(settelments_colormap); title('Original')


for epoch_idx = 1 : 1
    [fields_data_out] = adjust_frontiers_029(fields_data_in, basic_penetrating_element_shape, epoch_idx, show)
    figure(55 + epoch_idx-1); image(mod(fields_data_out.fields, size(settelments_colormap, 1)-4)); colormap(settelments_colormap); title(sprintf('Output at epoch %d ', epoch_idx));
  %  fields_data_in = fields_data_out
end;


