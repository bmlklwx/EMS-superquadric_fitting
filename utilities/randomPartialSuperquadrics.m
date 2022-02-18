function [point_partial, nomral_partial] = randomPartialSuperquadrics(x, arclength, percentage)

[point, normal_gt] = sphericalProduct_sampling(x, arclength);
num_pt = size(point, 2);
num_rand = floor(num_pt * percentage);
idx = randi(num_pt);
[mIdx, ~] = knnsearch(point', point(:, idx)', 'K', num_rand);
point_partial = point(:, mIdx);
nomral_partial = normal_gt(:, mIdx);
end