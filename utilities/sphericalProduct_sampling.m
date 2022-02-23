function [point] = sphericalProduct_sampling(para, arclength)
% Sampling of Superquadrics by spherical products (almost-uniform)
% 2D superellipse sampling referenced from
% Based on Pilu and Fisher, Equal-Distance Sampling of Superellipse Models, 1995

% Written and modified by Weixiao Liu, PhD student @ JHU, NUS
% June 2th, 2021, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION:

% INPUT:


% OUTPUT:

% -------------------------------------------------------------------------
% avoiding numerical instability when sampling points
epsilon = max(para(1 : 2), [0.007, 0.007]);

a = para(3 : 5);
R = eul2rotm(para(6 : 8));
t = para(9 : 11);

[~, theta1, critical1] = uniformSampledSuperellipse(epsilon(1), [(a(1) + a(2)) / 2, a(3)], [0 0 0], arclength, 0);

a_mod = zeros(1, size(theta1, 2));
a_mod(1 : critical1 - 1) = cos(theta1(1 : critical1 - 1)) .^ epsilon(1);
a_mod(critical1 : end) = flip(sin(theta1(critical1 : end))) .^ epsilon(1);

points_fw1 = angle2points(theta1(1 : critical1 - 1), [1, a(3)], epsilon(1));
points_bw1 = flip(angle2points(theta1(critical1 : end), [a(3), 1], epsilon(1)), 2);

point1 = [points_fw1, [points_bw1(2, :); points_bw1(1, :)]];

point2 = cell(1, size(a_mod, 2));

for i = 1 : size(a_mod, 2)
    [~, theta2, critical2] = uniformSampledSuperellipse(epsilon(2), [a(1), a(2)], [0 0 0], arclength / a_mod(i), 0);
    
    points_fw2 = angle2points(theta2(1 : critical2 - 1), [a(1), a(2)], epsilon(2));
    points_bw2 = flip(angle2points(theta2(critical2 : end), [a(2), a(1)], epsilon(2)), 2); 
    point2{i} = [points_fw2, [points_bw2(2, :); points_bw2(1, :)]];
end

% save in cell to differentiate points with their latitude
point_upper = cell(size(point1, 2), 1);

for i = 1 : size(point1, 2)
    point_temp = [point2{i} * point1(1, i); ones(1, size(point2{i}, 2)) * point1(2, i)];
    
    num_pt = size(point_temp, 2);
    point_upper{i} = [point_temp, ...
        [-point_temp(1, 1 : num_pt - 1); point_temp(2, 1 : num_pt - 1); point_temp(3, 1 : num_pt - 1)], ...
        [-point_temp(1, 2 : end); -point_temp(2, 2 : end); point_temp(3, 2 : end)], ...
        [point_temp(1, 2 : num_pt - 1); -point_temp(2, 2 : num_pt - 1); point_temp(3, 2 : num_pt - 1)]];   
end

point_upper{size(point1, 2)} = point_upper{size(point1, 2)}(:, 1);

num_pt = 0;

for i = 2 : size(point1, 2)
    point(:, num_pt + 1 : num_pt + size(point_upper{i}, 2)) = ...
        point_upper{i};
    num_pt = num_pt + size(point_upper{i}, 2);
end

point = [point_upper{1}, point, [point(1, :); point(2, :); -point(3, :)]]';

point = R * point' + t';
% -------------------------functions --------------------------------------

%------------------ mapping from angles to points -------------------------
    function [point] = angle2points(theta, scale, sigma)
        point = zeros(2, size(theta, 2));
        point(1, :) = scale(1) .* sign(cos(theta)) .* abs(cos(theta)).^sigma;
        point(2, :) = scale(2) .* sign(sin(theta)) .* abs(sin(theta)).^sigma;
    end

%------------------ sampling on superellipse uniformly --------------------
    function [point, theta, critical] = uniformSampledSuperellipse(sigma, scale, xform, arclength ,disp)
        
        threshold = 1e-2;
        num_limit = 10000;
        theta = zeros(1, num_limit);
        seg = zeros(1, num_limit);
        theta(1) = 0;
        seg(1) = 1;
        for m = 2 : num_limit
            [dt, seg_temp] = dtheta(theta(m - 1), arclength, threshold, scale, sigma);
            theta_temp = theta(m - 1) + dt;
            
            if theta_temp > pi/4
                break
            else
                if m < num_limit
                    theta(m) = theta_temp;
                    seg(m) = seg_temp;
                else
                    error(['The number of the sampled points exceeds the limit of ', ...
                        num2str(num_limit * 4),...
                        '. Please increase the arclength or raise the limit'])
                end
            end
        end
        critical = m;
        seg(critical) = 1;
        for n = critical + 1 : num_limit
            [dt, seg_temp] = dtheta(theta(n - 1), arclength, threshold, flip(scale), sigma);
            theta_temp = theta(n - 1) + dt;
            
            if theta_temp > pi/4
                break
            else
                if n < num_limit
                    theta(n) = theta_temp;
                    seg(n) = seg_temp;
                else
                    error(['The number of the sampled points exceeds the limit of ', ...
                        num2str(num_limit * 4),...
                        '. Please increase the arclength or raise the limit'])
                end
            end
        end
        
        num_point = n - 1;
        theta = theta(1 : num_point);
        seg = seg(1  : num_point);
        seg = [seg(1  : critical - 1), flip(seg(critical : end))];
        
        points_fw = angle2points(theta(1 : critical - 1), scale, sigma);
        points_bw = flip(angle2points(theta(critical : end), flip(scale), sigma), 2);
        point = [points_fw, [points_bw(2, :); points_bw(1, :)]];
        
        point = [point, [-point(1, 1 : num_point - 1); point(2, 1 : num_point - 1)], ...
            [-point(1, 2 : end); -point(2, 2 : end)], [point(1, 2 : num_point - 1); ...
            -point(2, 2 : num_point - 1)]];
        
        seg = [seg, seg(1 : num_point - 1), seg(2 : end), seg(2 : num_point - 1)];
        
        point = [cos(xform(1)), -sin(xform(1)); sin(xform(1)), cos(xform(1))]...
            * point + [xform(2); xform(3)];
        
        if disp == 1
            figure
            plot(point(1, seg == 1), point(2, seg == 1), '*')
            hold on
            plot(point(1, seg == 2), point(2, seg == 2), '*')
            hold off
            axis equal
        end
    end

%------------------ calculation of theta interval----- --------------------
    function [dt, seg] = dtheta(theta, arclength, threshold, scale, sigma)
        if theta < threshold
            dt = abs((arclength / scale(2) + (theta)^(sigma))^(1 / sigma) ...
                - (theta));
            seg = 1;
        else
            dt = arclength / sigma * ((cos(theta) ^ 2 * sin(theta) ^ 2) / ...
                (scale(1) ^ 2 * cos(theta) ^ (2 * sigma) * sin(theta) ^ 4 + ...
                scale(2) ^ 2 * sin(theta) ^ (2 * sigma) * cos(theta) ^ 4))^(1 / 2);
            seg = 2;
        end
    end

end