%
% compute smallest angle that could be obtained by chance for a given P,
% for empirically generated distributions of angles (computed with
% empirical_angle_distribution.m)
%
%   function angle_non_orth = get_random_angle_for_Pval( dist_angles, ...
%       space_dim, plane_dim, samples, Pval )
%
%

function angle_non_orth = get_random_angle_for_Pval( dist_angles, space_dim, ...
                            plane_dim, samples, Pval )


% turn the angles into a distribution of angles
hist_x              = 0:pi/2/90:pi/2;
hist_dist_angles    = zeros(length(hist_x)-1,length(plane_dim),length(space_dim));

for n = 1:length(space_dim)
    for d = 1:length(plane_dim)
        hist_dist_angles(:,d,n) = histcounts(dist_angles(:,d,n), hist_x)/samples;
    end
end

% find the largest angle with P > Pval
angle_non_orth      = zeros(length(plane_dim),length(space_dim));
for n = 1:length(space_dim)
    for d = 1:length(plane_dim)
        angle_non_orth(d,n) = find(cumsum(hist_dist_angles(:,d,n))>Pval, 1);
    end
end