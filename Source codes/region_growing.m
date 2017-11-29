function [fields_data, colored_fields] = region_growing(fields_data_in,weight)
% INOUT ===================================================================
% fields_data_in.im_d
% fields_data_in.im_file
% fields_data_in.places_need_to_be_cehcked_for_seeding
% fields_data_in.seeds_data{seeds_data_idx}.seeds_
% fields_data_in.seeds_data{seeds_data_idx}.seeds_{1}.equations_of_plan_
% fields_data_in.seeds_data{seeds_data_idx}.seeds_{1}.points
% fields_data_in.seeds_data{seeds_data_idx}.seeds_{1}.divergance_from_plan
% fields_data_in.seeds_data{seeds_data_idx}.initial_seeds_shape
% fields_data_in.seeds_data{seeds_data_idx}.min_allowed_points_for_initial_seeding
% fields_data_in.seeds_data{seeds_data_idx}.im_d
% fields_data_in.seeds_data{seeds_data_idx}.places_to_cehck_for_initial_seeding
% fields_data_in.seeds_data{seeds_data_idx}.divergance_from_plan_mat
% fields_data_in.seeds_data{seeds_data_idx}.order_of_generating_seeds
% fields_data_in.seeds_data{seeds_data_idx}.im_d = []; % just to save space; since the image will be already stored in "fields_data_in.im_d"

% fields_data_in.order_of_using_seeds

% OUTPUT ==================================================================

% fields_data (output) = fields_data_in (input) AND THE FOLLOWING FIELDS :

% fields_data_in.field_{1}.initial_seed.field_id = [];
% fields_data_in.field_{1}.initial_seed.positions = [];
% fields_data_in.field_{1}.points_of_filed = [];
% fields_data_in.field_{1}.equations_of_plan_{1} = [];
% fields_data_in.field_{1}.extent_of_grouth = [];
% fields_data_in.basic_data.maximum_allowed_dissimilarity
% fields_data_in.basic_data.lag_of_update_THR
% fields_data_in.idxs_of_using_seeds
% fields_data.seed_of_flat_area_mat: this figure shows all the seeds used to generate the flat areas, "the final used seeds figure"
%
% =========================================================================

global im_d;
global places_need_to_be_cehcked_for_seeding;
global mask_of_current_batch;
global fields;

fields_data = fields_data_in; 
im_d = fields_data.im_d;
[r,c] = size(im_d);
fields = 0*ones(r, c);
seeds_data_idx = length(fields_data.seeds_data);
fields_data.basic_data.current_idx_of_flat_area = 1;
fields_data.seed_of_flat_area_mat = 0*ones(r, c);
% Is it suitable for seeding
[places_need_to_be_cehcked_for_seeding] = where_to_start_seeding(fields_data); 
fields_data.places_need_to_be_cehcked_for_seeding = places_need_to_be_cehcked_for_seeding;
[fields_data] = get_order_of_using_seeds(fields_data, seeds_data_idx); 

for ii = 1: length(fields_data.idxs_of_using_seeds)
    ss = fields_data.idxs_of_using_seeds(ii);
    mask_of_current_batch = 0*ones(r, c); 
    points = fields_data.seeds_data{seeds_data_idx}.seeds_{ss}.points;
    idxs_of_initial_seed = sub2ind( size(places_need_to_be_cehcked_for_seeding), points(:, 1), points(:, 2) );
    valid_initial_seed = (sum( 0 == places_need_to_be_cehcked_for_seeding( idxs_of_initial_seed ))) == 0;
    
    if valid_initial_seed
        field_idx = fields_data.basic_data.current_idx_of_flat_area;
        bridth_of_grow = 0;
        places_need_to_be_cehcked_for_seeding(idxs_of_initial_seed) = 0;
        fields(idxs_of_initial_seed) = field_idx;
        mask_of_current_batch(idxs_of_initial_seed) = 1;
        fields_data.seeds_data{seeds_data_idx}.used_seeds(field_idx) = ss;
        fields_data.seed_of_flat_area_mat(idxs_of_initial_seed) = field_idx;
        fields_data.field_{field_idx}.seed = fields_data.seeds_data{seeds_data_idx}.seeds_{ss};
        fields_data.field_{field_idx}.points_of_filed = points;
        fields_data.field_{field_idx}.equations_of_plan_{1} = fields_data.seeds_data{seeds_data_idx}.seeds_{ss}.equations_of_plan_{1};
        set(0,'RecursionLimit',1500)  %for large surface growing
        [fields_data.field_{field_idx}.equations_of_plan_, new_points_of_filed, bridth_of_grow] =  grow_flat_area(fields_data.basic_data, fields_data.field_{field_idx}.equations_of_plan_, points, bridth_of_grow,weight,im_d);
        fields_data.field_{field_idx}.points_of_filed = vertcat(fields_data.field_{field_idx}.points_of_filed, new_points_of_filed);
        fields_data.field_{field_idx}.extent_of_grouth = bridth_of_grow;
        fields_data.basic_data.current_idx_of_flat_area = field_idx +1; % the final and real value of the current field is current_idx_of_flat_area-1
    end
end

fields_data.fields = fields;
fields_data.non_flat_area = (fields_data.fields == 0);
[r,c] = size(fields_data.fields);
colored_fields = zeros(r, c, 3);

if max(max(fields_data.fields)) < 0.95 * 2^24
    tmp_im = 0.95 * 2^24 * ( fields_data.fields / max(max(fields_data.fields)) );
    tmp_im = uint64(round(tmp_im -10));
    colored_fields(:, :, 1) = bitand(tmp_im,  2^8-1); %bit-wise AND; bit by bit AND
    colored_fields(:, :, 2) = bitand(bitshift(tmp_im, -8),  2^8-1);%(jz) bitshift(tmp_im, -8), the bits in tmo_im shift to right by 8 bits
    colored_fields(:, :, 3) = bitand(bitshift(tmp_im, -16),  2^8-1);
end

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 
function [fields_data] = get_order_of_using_seeds(fields_data, seeds_data_idx)
% this function is used to make sure the way of get initial seeds is the same way of using the seeds and also define how to realise the order

if (lower(fields_data.seeds_data{seeds_data_idx}.order_of_generating_seeds) ~= lower('RasterScan'))
    error('TTILLO: for now the only supported input is RasterScan for order_of_generating_seeds')
end

seeds_num = length(fields_data.seeds_data{seeds_data_idx}.seeds_);

switch lower(fields_data.order_of_using_seeds)
    case {'random'}
        fields_data.idxs_of_using_seeds = randperm(seeds_num);
        disp('Method is random')
        
    case 'rasterscan_col'
        fields_data.idxs_of_using_seeds = [1 : seeds_num];
        disp('Method is raster scan col')
        
    case 'rasterscan_row'
        tmp_idx_mat = fields_data.seeds_data{seeds_data_idx}.divergance_from_plan_mat;
        tmp_idx_mat( ~isnan(tmp_idx_mat) ) = [1 : seeds_num];
        tmp_idx_vct = tmp_idx_mat';
        tmp_idx_vct = tmp_idx_vct(:);
        fields_data.idxs_of_using_seeds = tmp_idx_vct(~isnan(tmp_idx_vct))';
        disp('Method is raster scan row')
        
    case 'ascendorder'  %raising
        divergance_from_plan_vct = zeros(1, seeds_num);
        for ii = 1 : seeds_num
            divergance_from_plan_vct(ii) = fields_data.seeds_data{seeds_data_idx}.seeds_{ii}.divergance_from_plan;
        end
        [~, Itmp] = sort(divergance_from_plan_vct, 'ascend');
        fields_data.idxs_of_using_seeds = Itmp;
        disp('Method is ascend order')
        
    case 'descendorder'
        divergance_from_plan_vct = zeros(1, seeds_num);
        for ii = 1 : seeds_num
            divergance_from_plan_vct(ii) = fields_data.seeds_data{seeds_data_idx}.seeds_{ii}.divergance_from_plan;
        end
        [~, Itmp] = sort(divergance_from_plan_vct, 'descend');
        fields_data.idxs_of_using_seeds = Itmp;
        disp('Method is descend order')
        
    case 'debug'
        divergance_from_plan_vct = zeros(1, seeds_num);
        for ii = 1 : seeds_num
            divergance_from_plan_vct(ii) = fields_data.seeds_data{seeds_data_idx}.seeds_{ii}.divergance_from_plan;
        end
        [~, Itmp] = sort(divergance_from_plan_vct, 'ascend');
        fields_data.idxs_of_using_seeds = Itmp;
        disp('Method is debug order')
        
        fields_data.idxs_of_using_seeds = Itmp(255: 555);
        
        
        
    otherwise
        error('TTILLO: non supported order of using seeds')
end

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



function [equations_of_plan_, points_of_filed_new, bridth_of_grow] = grow_flat_area(basic_data, equations_of_plan_, points_of_filed, bridth_of_grow,weight,im_d)
%this function is used to
global places_need_to_be_cehcked_for_seeding;
global mask_of_current_batch;
global fields;

field_idx = basic_data.current_idx_of_flat_area;
points_of_filed_new = points_of_filed;
[neighbors_points, valid_neighbors] = get_neighbors_points();%(jz)here there is no any judgement on the neighbor points expect whether it is zero or not.
if ~valid_neighbors
    return;
else
    
    [D_square_vct, points_belong_to_surfnace, similarity_THR] = get_points_distances_from_plan(basic_data, equations_of_plan_, neighbors_points, bridth_of_grow,weight);    
    [I_vct] = find(points_belong_to_surfnace == 1); % I_vct is the index of colume scan results
    if ~isempty(I_vct) % i_vct is not empty then ...
        bridth_of_grow = bridth_of_grow +1;
        % in neighbors_points the values of first col and second col will be chosen based on the value of I_vct
        indexis_tmp = sub2ind(size(fields), neighbors_points(I_vct, 1), neighbors_points(I_vct, 2)); % these positions will find the index on "field"
        fields(indexis_tmp) = field_idx;
        places_need_to_be_cehcked_for_seeding(indexis_tmp) = 0;
        mask_of_current_batch(indexis_tmp) = 1;
        
        points_of_filed_new = vertcat(points_of_filed, neighbors_points(I_vct, :)); %(jz) vertically concatenate two matrix
        [equations_of_plan_] = refine_plan_parameters(basic_data, equations_of_plan_, points_of_filed_new, bridth_of_grow);
        [equations_of_plan_, points_of_filed_new, bridth_of_grow] = grow_flat_area(basic_data, equations_of_plan_, points_of_filed_new, bridth_of_grow,weight,im_d);
    else
        return;
    end
end
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [neighbors_points, valid_neighbors] = get_neighbors_points()

global im_d;
global places_need_to_be_cehcked_for_seeding;
global mask_of_current_batch;


mask_neighbors = places_need_to_be_cehcked_for_seeding & xor(mask_of_current_batch, imdilate(mask_of_current_batch, ones(3)));

figure(4); image(uint8(mask_neighbors)); colormap(gray(2));
[ri_idx, ci_idx] = find(mask_neighbors == 1);

neighbors_points( : , 1 ) = ri_idx;
neighbors_points( : , 2 ) = ci_idx;
neighbors_points( : , 3 ) = im_d(mask_neighbors == 1);

valid_neighbors =  ~isempty(ri_idx);

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



function [places_need_to_be_cehcked_for_seeding_final] = where_to_start_seeding(fields_data)
% this function is used to find a nonzero area in im_d as the seed to start
places_need_to_be_cehcked_for_seeding_final = (fields_data.im_d ~= 0);
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [equations_of_plan_] = refine_plan_parameters(basic_data, equations_of_plan_, points_of_filed, bridth_of_grow)

field_idx = basic_data.current_idx_of_flat_area;
[A] = equations_of_plan_{ length(equations_of_plan_) }(1);
[B] = equations_of_plan_{ length(equations_of_plan_) }(2);
[C] = equations_of_plan_{ length(equations_of_plan_) }(3);
[D] = equations_of_plan_{ length(equations_of_plan_) }(4);


ri_idx = points_of_filed(:, 1);
ci_idx = points_of_filed(:, 2);
z_data = points_of_filed(:, 3);
sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare');

% [ A, B, C, D] = [ sf.p10, sf.p01, -1, sf.p00];
%     [ sf.p10, sf.p01, -1, sf.p00]./[ A, B, C, D]
%     [ A, B, C, D]
d_new = sum((z_data - feval(sf, ri_idx, ci_idx)).^2); %(jz) "feval(sf, ri_idx, ci_idx)" evaluate the function in sf when variables are ri_idx and ci_idx

sf_or = sf;
sf_or.p10 = A/(-C);
sf_or.p01 = B/(-C);
sf_or.p00 = D/(-C);
d_old = sum((z_data - feval(sf_or, ri_idx, ci_idx)).^2);

A =  sf.p10;
B = sf.p01;
C = -1;
D = sf.p00;
equations_of_plan_{ length(equations_of_plan_) +1 } = [A, B, C, D, bridth_of_grow];

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


function [D_square_vct, points_belong_to_surfnace,similarity_THR] = get_points_distances_from_plan(basic_data, equations_of_plan_, points, bridth_of_grow,weight)

[A] = equations_of_plan_{ length(equations_of_plan_) }(1);
[B] = equations_of_plan_{ length(equations_of_plan_) }(2);
[C] = equations_of_plan_{ length(equations_of_plan_) }(3);
[D] = equations_of_plan_{ length(equations_of_plan_) }(4);

D_square_vct = (points * [A B C]' + D).^2 / ([A B C] *[A B C]');

[similarity_THR] = get_learning_weight(basic_data, bridth_of_grow +1,points,weight);

points_belong_to_surfnace = (D_square_vct <= similarity_THR);

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

function [similarity_THR] = get_learning_weight(basic_data, bridth_of_grow,points,weight)
global im_d;
learning_w_per_point = 1 - exp(-bridth_of_grow./basic_data.lag_of_update_THR);
if bridth_of_grow < floor(size(im_d,1)*size(im_d,2)/400)
    similarity_THR = (basic_data.maximum_allowed_dissimilarity * learning_w_per_point).^2  ;
else
    similarity_THR = (weight*(255-points(:,3)).*basic_data.maximum_allowed_dissimilarity*learning_w_per_point).^2;
end
return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


