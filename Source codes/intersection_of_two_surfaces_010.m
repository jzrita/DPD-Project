function [plans_intersection, mask_of_intersection] = intersection_of_two_surfaces_010(fields_data, idx_of_plan1, idx_of_plan2, show, shift_from_boarder)

show_im = show.show_im;
fig_idx = show.fig_idx;
settelments_colormap = show.settelments_colormap;

equations_of_plan1 = fields_data.field_{idx_of_plan1}.equations_of_plan_{length(fields_data.field_{idx_of_plan1}.equations_of_plan_)};
equations_of_plan2 = fields_data.field_{idx_of_plan2}.equations_of_plan_{length(fields_data.field_{idx_of_plan2}.equations_of_plan_)};

[A] = equations_of_plan1(1);
[B] = equations_of_plan1(2);
[C] = equations_of_plan1(3);
[D] = equations_of_plan1(4);
n1 = [A B C];
d1 = -D;

[A] = equations_of_plan2(1);
[B] = equations_of_plan2(2);
[C] = equations_of_plan2(3);
[D] = equations_of_plan2(4);
n2 = [A B C];
d2 = -D;

% points_cordinates = pp + rr* tt; The equation of the intersection line between two surfaces
rr = cross(n1, n2); %gives the normal vector of the intersection line
n1_dot_n2 = n1 * n2';
denum = n1*n1' * n2*n2' - (n1_dot_n2)^2;
pp = ((d1*n2*n2' - d2* n1_dot_n2)*n1 + (d2*n1*n1' - d1* n1_dot_n2)*n2) / denum; 

%this above equation has exact same results as the below one(20150129)
% upper1= (d1*n2-d2*n1);
% upper2=cross(n1,n2);
% upper=cross(upper1,upper2);
% pp0=upper/(abs(cross(n1,n2))*abs(cross(n1,n2))')

% finding the limits 
tt_limit1 = ((1 +shift_from_boarder) - pp(1))/rr(1);
tt_limit2 = ((size(fields_data.fields, 1) -shift_from_boarder) - pp(1))/rr(1);
tt_min1 = min(tt_limit1, tt_limit2);
tt_max1 = max(tt_limit1, tt_limit2);

tt_limit1 = ((1 +shift_from_boarder) - pp(2))/rr(2);
tt_limit2 = ((size(fields_data.fields, 2) -shift_from_boarder) - pp(2))/rr(2);
tt_min2 = min(tt_limit1, tt_limit2);
tt_max2 = max(tt_limit1, tt_limit2);

im = fields_data.fields;
mask_of_intersection = zeros(size(im));
if (tt_max2 <= tt_min1) | (tt_max1 <= tt_min2)
    plans_intersection.valid_intersection = 0; 
else
    plans_intersection.valid_intersection = 1;
    tt_min = max(tt_min1, tt_min2);
    tt_max = min(tt_max1, tt_max2);
    
    tt = linspace(tt_min, tt_max, 4*max(size(fields_data.fields)))';
    points = [(pp(1) + tt*rr(1))  ,  (pp(2) + tt*rr(2))  ,  (pp(3) + tt*rr(3))];
    if ~isnan(points)
    points = round(points);
    [points, unique_idx]  = unique(points, 'rows');
    tt = tt(unique_idx);
    im(sub2ind(size(im), points(:, 1), points(:, 2))) = 30;
    mask_of_intersection(sub2ind(size(im), points(:, 1), points(:, 2))) = 1;
    plans_intersection.points = points;
    
    % points = pp + rr* t; The equation of the intersection line between two surfaces
    plans_intersection.tri_d_line_equation.pp = pp;
    plans_intersection.tri_d_line_equation.rr = rr;
    plans_intersection.tri_d_line_equation.tt = tt;
    
    % points = pp + rr* x; The equation of the intersection line (in 2D) between two surfaces
    point_tmp = pp + rr*(0 - pp(1))/rr(1);
    plans_intersection.two_d_line_equation.pp = point_tmp(2);
    plans_intersection.two_d_line_equation.rr = rr(2)/rr(1); 
    else
        mask_of_intersection=nan;
    end
end;

% if show_im
%     figure(fig_idx + 5); image(fields_data.fields); colormap(settelments_colormap); title('Original');
%     if plans_intersection.valid_intersection
%         figure(fig_idx + 6); image(im); colormap(settelments_colormap); title(sprintf('With the intersection line between two surfaces'));
%     else
%         figure(fig_idx + 6); image(im); colormap(settelments_colormap); title(sprintf('NO intersection between the two surfaces in the viewed image'));
%     end;
% end

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


show.show_im = 1;
show.fig_idx = 30;
load('settelments_colormap.mat')

idx_of_plan1 = 6;
idx_of_plan2 = 28;

shift_from_boarder = 0;

[plans_intersection, mask_of_intersection] = intersection_of_two_surfaces_010(fields_data, idx_of_plan1, idx_of_plan2, show, shift_from_boarder);

% 
% for ii = 1 : numel(fields_data.field_)