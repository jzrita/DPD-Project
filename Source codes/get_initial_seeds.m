%% GO through (run) this code step by step and see how each parameter change.20131231


function  seeds_data = get_initial_seeds(seeds_data_in, im_base_name)
%%----------The output------------------%%%
% seeds_data.seeds_
% seeds_data.initial_seeds_shape
% seeds_data.min_allowed_points_for_initial_seeding
% seeds_data.im_d
% seeds_data.places_to_cehck_for_initial_seeding
% seeds_data.divergance_from_plan_mat

im_d = seeds_data_in.im_d;
[r, c] = size(im_d); % the size of the input image

places_to_cehck_for_initial_seeding = seeds_data_in.places_to_cehck_for_initial_seeding;
initial_seeds_shape = seeds_data_in.initial_seeds_shape;
min_allowed_points_for_initial_seeding = seeds_data_in.min_allowed_points_for_initial_seeding; %(jz)it's decided by the seeds shape
seed_shape = seeds_data_in.seed_shape;
[r_seed, c_seed] = size(initial_seeds_shape);
divergance_from_plan_mat = NaN * ones(r, c); %(initialise the matrix of divergance plan with value NaN)
idx = 0;
[ri_idx_base, ci_idx_base] = find( initial_seeds_shape ~= 0 );

for ci = 1 : (c - c_seed + 1) % don't change the order of scaning the matrix %(jz) the scaning order is from top to bottom and then from left to right
    fprintf('Processing column %d\n', ci);
    for ri = 1 : (r - r_seed + 1) %candidate_points can be treated as a butter
        candidate_points( : , 1 ) = ri_idx_base +ri -1;
        candidate_points( : , 2 ) = ci_idx_base +ci -1;
        %using the function 'sub2ind' to find the exact point in the im_d
        candidate_points( : , 3 ) = im_d( sub2ind(size(im_d), candidate_points( : , 1 ) , candidate_points( : , 2 )) ); % if the input image is depth map then it's Z value
        [valid_points, candidate_points] = are_points_valid_candidate_for_initial_seeding(candidate_points, places_to_cehck_for_initial_seeding, min_allowed_points_for_initial_seeding);
        
        if valid_points
            idx = idx +1;
            %each valid point will become a seed and has its own equations_of_plan
            [seeds_data.seeds_{idx}.equations_of_plan_, distances_from_plan_vct] = get_best_fitting_plan_for_points(candidate_points);
            seeds_data.seeds_{idx}.points = candidate_points;
            divergance_from_plan_mat(ri, ci) = mean(distances_from_plan_vct); %(jz)find the average divergance of each point in the seeds shape to the fitting surface.
            seeds_data.seeds_{idx}.divergance_from_plan = divergance_from_plan_mat(ri, ci);
        end
        clear candidate_points;
    end
end

seeds_data.initial_seeds_shape = seeds_data_in.initial_seeds_shape;
seeds_data.min_allowed_points_for_initial_seeding = seeds_data_in.min_allowed_points_for_initial_seeding;
seeds_data.im_d = im_d;
seeds_data.places_to_cehck_for_initial_seeding = places_to_cehck_for_initial_seeding;
seeds_data.divergance_from_plan_mat = divergance_from_plan_mat;
seeds_data.order_of_generating_seeds = seeds_data_in.order_of_generating_seeds;
seeds_data.seed_shape = seed_shape;
show_im=1;
fig_idx=1;
if show_im
    figure(fig_idx); image(double(divergance_from_plan_mat)); colormap(gray(255)); title(sprintf('the divergance from plan for seeds'));
end

    results_file = sprintf('%s_SeedName_%s_SeedSize_%03d.mat', im_base_name, upper(seeds_data.seed_shape), sum(sum(initial_seeds_shape)));
    save(results_file);

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

%%
function [valid_points, candidate_points_2nd_filtered] = are_points_valid_candidate_for_initial_seeding(candidate_points, places_to_cehck_for_initial_seeding, min_allowed_points_for_initial_seeding)

idx_valid_points_1st_filtered = (candidate_points(:, 3) ~= 0); %all the index of the nonzero point intensity in "idx_valid_points_1st_filtered"
candidate_points_1st_filtered =  candidate_points(idx_valid_points_1st_filtered, :); %all the onozero point intensity in candidate_points buffer, e.g. candidate_points (1,:)= the first row of candidate_points
idx_valid_points_2nd_filtered  = find( places_to_cehck_for_initial_seeding(  sub2ind( size(places_to_cehck_for_initial_seeding), candidate_points_1st_filtered(:, 1), candidate_points_1st_filtered(:, 2) ) ) );
candidate_points_2nd_filtered =  candidate_points_1st_filtered(idx_valid_points_2nd_filtered, :); %candidate_points_1st_filtered (1,:)= the first row of candidate_points_1st_filtered
valid_points = ( size(candidate_points_2nd_filtered, 1) >= min_allowed_points_for_initial_seeding); 

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||



%%
function [equations_of_plan_, distances_from_plan_vct] = get_best_fitting_plan_for_points(candidate_points)

%    d = -ax_0-by_0-cz_0
ri_idx = candidate_points(:, 1);
ci_idx = candidate_points(:, 2);
z_data = candidate_points(:, 3);
sf = fit( [ri_idx, ci_idx], z_data, 'poly11', 'Robust', 'Bisquare'); %Matlab built-in function
% distances_from_plan_vct is the distance vector of each point to the fitting surface
distances_from_plan_vct = (z_data - feval(sf, ri_idx, ci_idx)).^2; %feval is Matlab built-in function;

%      Linear model Poly11:
%      sf(x,y) = p00 + p10*x + p01*y
%      z = p00 + p10*x + p01*y
%     (jz) 0=  p10*x + p01*y-z+ p00
%      p00 = - p10*x - p01*y + z
%   [ A, B, C, D] = [ sf.p10, sf.p01, -1, sf.p00];

A =  sf.p10;
B = sf.p01;
C = -1;   %(jz) the coefficient of z
D = sf.p00;
equations_of_plan_{ 1 } = [A, B, C, D, 0]; %(jz) the "0" here is left for future recording the update times, since it could have have more fitting surface functions.

return %||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||


