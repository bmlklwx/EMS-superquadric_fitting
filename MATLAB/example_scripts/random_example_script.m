%% Generate superquadric and sampling points
clear
close all

% generate a random superquadric point cloud
% generating superquadric parameters
epsilon = 2 * rand(1, 2) + 0.004;
a = 2 * rand(1, 3) + 0.5;
euler = [2 * pi * rand, 2 * pi * rand, pi * rand] ;
t = 0.2 * rand(1, 3) - 0.1;
R = eul2rotm(euler);
x_gt = [epsilon, a, euler, t];

% point cloud sampling arclength
arclength = 0.2;

% generating points and the corresponding ground truth surface normal direction 
[point] = sphericalProduct_sampling(x_gt, arclength);
figure(1)
showPoints(point)
title('Original Point Cloud')

% generating random partial view
partial_ratio = 0.6;
[point] = randomPartialSuperquadrics(x_gt, arclength, partial_ratio);
figure(2)
showPoints(point)
hold on
title('Partial Point Cloud')

% add noise
noise_level = 0.1;
noise = rand(3, size(point, 2)) * a(1) * noise_level - a(1) * noise_level / 2;
point = point + noise;
figure(3)
showPoints(point)
axis equal
title('Partial Point Cloud with Noise')

% add outlier
outlier_ratio = 0.4;
outlier = mvnrnd(t, 2 .* eye(3), floor(outlier_ratio * size(point, 2)))';
point = [point, outlier];

figure(4)
showPoints(point)
axis equal
title('Partial Point Cloud with Noise and Outliers')

%% Superquadric Recovery

[x_ns] = numerical_fitting(point);
[x_radial] = superquadricsFitting(point, 'Radial');
[x_ems] = EMS(point, 'OutlierRatio', 0.2);

disp('---------------------------------------------------------------------')
disp('Groud Truth parameters are');
disp(x_gt)
disp('---------------------------------------------------------------------')
disp('NS Fitted parameters are')
disp(x_ns)
disp('Radial-LSQ Fitted parameters are')
disp(x_radial)
disp('EMS Fitted parameters are')
disp(x_ems)
disp('---------------------------------------------------------------------')

figure(5)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_ns, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.1, 'Light', 1);
hold off
title('NS Recovered Superquadric')

figure(6)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_radial, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.1, 'Light', 1);
hold off
title('Radial-LSQ Recovered Superquadric')

figure(7)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_ems, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.1, 'Light', 1);
hold off
title('EMS Recovered Superquadric')