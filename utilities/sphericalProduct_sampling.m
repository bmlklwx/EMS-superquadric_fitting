function [point, normal, etaCorrPoint, omegaCorrPoint] = sphericalProduct_sampling(para, arclength)
% Sampling of Superquadrics by spherical products (like-uniform)
% 2D superellipse sampling referenced from
% Algorithm by Pilu and Fisher, Equal-Distance Sampling of Superellipse Models, 1995

% Written by Weixiao Liu, PhD student @ JHU, NUS
% June 2th, 2021, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION:

% INPUT:


% OUTPUT:

% -------------------------------------------------------------------------
epsilon = max(para(1 : 2), [0.007, 0.007]); %0.01
a = para(3 : 5);
R = eul2rotm(para(6 : 8));
t = para(9 : 11);

[~, theta1, critical1] = uniformSampledSuperellipse(epsilon(1), [(a(1) + a(2)) / 2, a(3)], [0 0 0], arclength, 0);

eta_flipped = [theta1(1 : critical1 - 1), flip(theta1(critical1 : end))];

a_mod = zeros(1, size(theta1, 2));
a_mod(1 : critical1 - 1) = cos(theta1(1 : critical1 - 1)) .^ epsilon(1);
a_mod(critical1 : end) = flip(sin(theta1(critical1 : end))) .^ epsilon(1);

points_fw1 = angle2points(theta1(1 : critical1 - 1), [1, a(3)], epsilon(1));
normals_fw1 = angle2normals(theta1(1 : critical1 - 1), a, epsilon(1));

points_bw1 = flip(angle2points(theta1(critical1 : end), [a(3), 1], epsilon(1)), 2);
normals_bw1 = flip(angle2normals(theta1(critical1 : end), [a(3), 1], epsilon(1)), 2);

point1 = [points_fw1, [points_bw1(2, :); points_bw1(1, :)]];
normal1 = [normals_fw1, [normals_bw1(2, :); normals_bw1(1, :)]];

point2 = cell(1, size(a_mod, 2));
normal2 = cell(1, size(a_mod, 2));
omega_flipped = cell(1, size(a_mod, 2));

for i = 1 : size(a_mod, 2)
    [~, theta2, critical2] = uniformSampledSuperellipse(epsilon(2), [a(1), a(2)], [0 0 0], arclength / a_mod(i), 0);
    omega_flipped{i} = [theta2(1 : critical2 - 1), flip(theta2(critical2 : end))];
    
    points_fw2 = angle2points(theta2(1 : critical2 - 1), [a(1), a(2)], epsilon(2));
    normals_fw2 = angle2normals(theta2(1 : critical2 - 1), [a(1), a(2)], epsilon(2));
    
    points_bw2 = flip(angle2points(theta2(critical2 : end), [a(2), a(1)], epsilon(2)), 2);
    normals_bw2 = flip(angle2normals(theta2(critical2 : end), [a(2), a(1)], epsilon(2)), 2);
    
    point2{i} = [points_fw2, [points_bw2(2, :); points_bw2(1, :)]];
    normal2{i} = [normals_fw2, [normals_bw2(2, :); normals_bw2(1, :)]];
end

% save in cell to differentiate points with their latitude
point_upper = cell(size(point1, 2), 1);
normal_upper = cell(size(point1, 2), 1);
eta_upper = cell(size(point1, 2), 1);
omega_upper = cell(size(point1, 2), 1);

for i = 1 : size(point1, 2)
    point_temp = [point2{i} * point1(1, i); ones(1, size(point2{i}, 2)) * point1(2, i)];
    normal_temp = [normal2{i} * normal1(1, i); ones(1, size(normal2{i}, 2)) * normal1(2, i)];
    
    eta_temp = eta_flipped(i) * ones(1, size(point2{i}, 2));
    omega_temp = omega_flipped{i};
    
    num_pt = size(point_temp, 2);
    point_upper{i} = [point_temp, ...
        [-point_temp(1, 1 : num_pt - 1); point_temp(2, 1 : num_pt - 1); point_temp(3, 1 : num_pt - 1)], ...
        [-point_temp(1, 2 : end); -point_temp(2, 2 : end); point_temp(3, 2 : end)], ...
        [point_temp(1, 2 : num_pt - 1); -point_temp(2, 2 : num_pt - 1); point_temp(3, 2 : num_pt - 1)]];
    normal_upper{i} = [normal_temp, ...
        [-normal_temp(1, 1 : num_pt - 1); normal_temp(2, 1 : num_pt - 1); normal_temp(3, 1 : num_pt - 1)], ...
        [-normal_temp(1, 2 : end); -normal_temp(2, 2 : end); normal_temp(3, 2 : end)], ...
        [normal_temp(1, 2 : num_pt - 1); -normal_temp(2, 2 : num_pt - 1); normal_temp(3, 2 : num_pt - 1)]];
    
    eta_upper{i} = [eta_temp, eta_temp(1 : end - 1), eta_temp(2 : end), eta_temp(2 : end - 1)];
    omega_upper{i} = [omega_temp, omega_temp(1 : end - 1), omega_temp(2 : end), omega_temp(2 : end - 1)];
    
end

point_upper{size(point1, 2)} = point_upper{size(point1, 2)}(:, 1);
normal_upper{size(point1, 2)} = normal_upper{size(point1, 2)}(:, 1);
eta_upper{size(point1, 2)} = 0;
omega_upper{size(point1, 2)} = 0;

num_pt = 0;

for i = 2 : size(point1, 2)
    point(:, num_pt + 1 : num_pt + size(point_upper{i}, 2)) = ...
        point_upper{i};
    normal(:, num_pt + 1 : num_pt + size(point_upper{i}, 2)) = ...
        normal_upper{i};
    etaCorrPoint(:, num_pt + 1 : num_pt + size(point_upper{i}, 2)) = eta_upper{i};
    omegaCorrPoint(:, num_pt + 1 : num_pt + size(point_upper{i}, 2)) = omega_upper{i};
    num_pt = num_pt + size(point_upper{i}, 2);
end

point = [point_upper{1}, point, [point(1, :); point(2, :); -point(3, :)]]';
etaCorrPoint = [eta_upper{1}, etaCorrPoint, etaCorrPoint]';
omegaCorrPoint = [omega_upper{1}, omegaCorrPoint, omegaCorrPoint]';
normal = [normal_upper{1}, normal, [normal(1, :); normal(2, :); -normal(3, :)]]';

point = R * point' + t';
normal = R * normal';

end