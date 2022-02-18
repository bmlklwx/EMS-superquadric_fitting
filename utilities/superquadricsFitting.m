function [x] = superquadricsFitting(point, cost_function_type)
%
%

% Written by Weixiao Liu, PhD student @ JHU, NUS
% May 24th, 2021, Singapore
% -------------------------------------------------------------------------
% DESCRIPTION:

% INPUT:


% OUTPUT:

% -------------------------------------------------------------------------

options = optimoptions('lsqnonlin', 'Algorithm', 'levenberg-marquardt', 'Display','off', 'MaxIterations', 50);
% 'StepTolerance', 1e-10
% 'trust-region-reflective','levenberg-marquardt', 'iter'

%% Initialization

t0 = mean(point, 2);
point = point - t0;
[EigenVector, ~] = EigenAnalysis(point);

EigenVector0 = [EigenVector(:, 1), EigenVector(:, 2), cross(EigenVector(:, 1), EigenVector(:, 2))];
EigenVector1 = [EigenVector(:, 1), EigenVector(:, 3), cross(EigenVector(:, 1), EigenVector(:, 3))];
EigenVector2 = [EigenVector(:, 2), EigenVector(:, 3), cross(EigenVector(:, 2), EigenVector(:, 3))];  
euler0 = [rotm2eul(EigenVector0); rotm2eul(EigenVector1); rotm2eul(EigenVector2)];

point_rot0 = EigenVector0' * point;
point_rot1 = EigenVector1' * point;
point_rot2 = EigenVector2' * point;

s0 = [median(abs(point_rot0(1, :))), median(abs(point_rot0(2, :))), median(abs(point_rot0(3, :)));
      median(abs(point_rot1(1, :))), median(abs(point_rot1(2, :))), median(abs(point_rot1(3, :)));
      median(abs(point_rot2(1, :))), median(abs(point_rot2(2, :))), median(abs(point_rot2(3, :)));];      
      
x0 = [ones(3, 2), s0, euler0, zeros(3, 3)];

upper = 4 * max(max(abs(point)));

lb = [0.1 0.1 0.01 0.01 0.01 -2*pi -2*pi -2*pi -ones(1, 3) * upper]; %0.05
ub = [2.0 2.0 ones(1, 3) * upper  2*pi 2*pi 2*pi ones(1, 3) * upper];

%% Function selection and Fitting

if strcmp(cost_function_type,'Vanilla_IO')
    cost_func = @(x) vanilla_in_out_func(x, point);
end

if strcmp(cost_function_type, 'Power_IO')
    cost_func = @(x) power_in_out_func(x, point);
end

if strcmp(cost_function_type, 'Vol_Power_IO')
    cost_func = @(x) volume_power_in_out_func(x, point);
end

if strcmp(cost_function_type, 'Radial')
    cost_func = @(x) radial_distance_func(x, point);
end

if strcmp(cost_function_type, 'Reproject')
    point = point(:, point(2, :) ~= 0);
    cost_func = @(x) reproject_distance_func(x, point);
end

if strcmp(cost_function_type, 'Normal')
    point = point(:, point(2, :) ~= 0);
    cost_func = @(x) normal_func(x, point);
end

x = zeros(3, 11);
residue = Inf * ones(1, 3);
for i = 1 : size(x0, 1)
    [x(i, :), residue(i)] = lsqnonlin(cost_func, x0(i, :),lb,ub,options);
end
[~, idx] = min(residue);
x = x(idx, :);

% rotation = EigenVector * eul2rotm(x(6 : 8));
% x(6 : 8) = rotm2eul(rotation);
% x(9 : 11) = (EigenVector * x(9 : 11)' + t0)';

x(9 : 11) = x(9 : 11) + t0';

% -------------------------------------------------------------------------
    function [value] = vanilla_in_out_func(para, X)

        R = eul2rotm(para(6 : 8));
        t = para(9 : 11);      
        
        X_c = R' * X - R' * t';
        %         X_c = X;
        value = ((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1)) - 1)';
    end

    function [value] = power_in_out_func(para, X)
        % Gross and Boult, 1988
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
         
        X_c = R' * X - R' * t';
        value = (((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ para(1) - 1)';
    end

    function [value] = volume_power_in_out_func(para, X)
        % Solina and Bajcsy, 1990
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
        
        X_c = R' * X - R' * t';
        value = (para(3) * para(4) * para(5)) ^ (1/ 2) .* (((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ para(1) - 1)';
    end

    function [value] = radial_distance_func(para, X)
        % Gross and Boult, 1988
        
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
        
        X_c = R' * X - R' * t';
        r_0 = vecnorm(X_c);
        dist = abs(((((X_c(1, :) / para(3)) .^ (2)) .^ (1 / para(2)) + ...
            ((X_c(2, :) / para(4)) .^ (2)) .^ (1 / para(2))) .^ (para(2) / para(1)) + ...
            ((X_c(3, :) / para(5)) .^ (2)) .^ (1 / para(1))) .^ (-para(1) / 2) - 1);
        value = r_0 .* dist;
        
    end

    function [value] = reproject_distance_func(para, X)
               
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
        
        X_c = R' * X - R' * t';
        [eta, omega, flag] = point2angle(X_c, para(1 : 2), para(3 : 5));
        [X_reproject] = angluarReprojection(eta, omega, para(1 : 2), para(3 : 5), flag);
        
        value = vecnorm(X_c - X_reproject);
    end
    function [value] = normal_func(para, X)
               
        R = eul2rotm(para(6 : 8));
        t = para(9 : 11); 
        
        X_c = R' * X - R' * t';
        [eta, omega, flag] = point2angle(X_c, para(1 : 2), para(3 : 5));
        
        [X_reproject] = angluarReprojection(eta, omega, para(1 : 2), para(3 : 5), flag);
        [normal_reproject] = angleflag2normal(eta, omega, para(1 : 2), para(3 : 5), flag);
%         value = abs(dot(X_c - X_reproject, normal_reproject));
        value = abs(dot(X_c - X_reproject, normal_reproject)) + vecnorm(X_c - X_reproject);
    end
end