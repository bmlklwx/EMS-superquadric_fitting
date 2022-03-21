%% loading point cloud
clear
close all

% add path
addpath('../src/')
addpath('../src/utilities/')

% read .ply
pc = pcread('./data/noisy_pointCloud_example_1.ply');
point = pc.Location';

%% Superquadric Recovery

% applying ns method
[x_ns] = numerical_fitting(point);
% applying radial-lsq method
[x_radial] = superquadricsFitting(point, 'Radial');
% applying EMS method
[x_ems] = EMS(point);

% print ground truth and recovered parameters
disp('---------------------------------------------------------------------')
disp('NS Fitted parameters are')
disp(x_ns)
disp('Radial-LSQ Fitted parameters are')
disp(x_radial)
disp('EMS Fitted parameters are')
disp(x_ems)
disp('---------------------------------------------------------------------')

% plot input points and rendering recovered superquadric
figure(1)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_ns, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.05, 'Light', 1);
hold off
title('NS Recovered Superquadric')

figure(2)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_radial, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.05, 'Light', 1);
hold off
title('Radial-LSQ Recovered Superquadric')

figure(3)
showPoints(point, 'Color', 'r')
hold on
showSuperquadrics(x_ems, 'Color', [0 0 1], 'FaceAlpha', 0.7, 'Arclength', 0.05, 'Light', 1);
hold off
title('EMS Recovered Superquadric')