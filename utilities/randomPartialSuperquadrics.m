function [point_partial] = randomPartialSuperquadrics(x, arclength, percentage)

[point] = sphericalProduct_sampling(x, arclength);
num_pt = size(point, 2);
num_rand = floor(num_pt * percentage);
idx = randi(num_pt);
[mIdx, ~] = knnsearch(point', point(:, idx)', 'K', num_rand);
point_partial = point(:, mIdx);

end