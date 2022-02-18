function [point, theta, critical] = uniformSampledSuperellipse(sigma, scale, xform, arclength ,disp)
% Uniform sampling of Superellipses
% Algorithm by Pilu and Fisher, Equal-Distance Sampling of Superellipse Models, 1995

% Written by Weixiao Liu, PhD student @ JHU, NUS
% May 24th, 2021, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION:

% INPUT:


% OUTPUT:

% -------------------------------------------------------------------------
threshold = 1e-2;
num_limit = 10000;
theta = zeros(1, num_limit);
seg = zeros(1, num_limit);
theta(1) = 0;
seg(1) = 1;
for i = 2 : num_limit
    [dt, seg_temp] = dtheta(theta(i - 1), arclength, threshold, scale, sigma);
    theta_temp = theta(i - 1) + dt;
    
    if theta_temp > pi/4
        break
    else
        if i < num_limit
            theta(i) = theta_temp;
            seg(i) = seg_temp;
        else
            error(['The number of the sampled points exceeds the limit of ', ...
                num2str(num_limit * 4),...
                '. Please increase the arclength or raise the limit'])
        end
    end
end
critical = i;
seg(critical) = 1;
for j = critical + 1 : num_limit
    [dt, seg_temp] = dtheta(theta(j - 1), arclength, threshold, flip(scale), sigma);
    theta_temp = theta(j - 1) + dt;
    
    if theta_temp > pi/4
        break
    else
        if j < num_limit
            theta(j) = theta_temp;
            seg(j) = seg_temp;
        else
            error(['The number of the sampled points exceeds the limit of ', ...
                num2str(num_limit * 4),...
                '. Please increase the arclength or raise the limit'])
        end
    end
end

num_pt = j;
theta = theta(1 : num_pt);
seg = seg(1  : num_pt);
seg = [seg(1  : critical - 1), flip(seg(critical : end))];

points_fw = angle2points(theta(1 : critical - 1), scale, sigma);
points_bw = flip(angle2points(theta(critical : end), flip(scale), sigma), 2);
point = [points_fw, [points_bw(2, :); points_bw(1, :)]];

point = [point, [-point(1, 1 : num_pt - 1); point(2, 1 : num_pt - 1)], ...
    [-point(1, 2 : end); -point(2, 2 : end)], [point(1, 2 : num_pt - 1); ...
    -point(2, 2 : num_pt - 1)]];

seg = [seg, seg(1 : num_pt - 1), seg(2 : end), seg(2 : num_pt - 1)];

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
