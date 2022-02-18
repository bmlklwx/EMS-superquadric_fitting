%% Generate superquadric and sampling points
clear
close all

addpath('/home/saintsbury/src/Superquadrics_Fitting/')
addpath('/home/saintsbury/src/Superquadrics_Fitting/functions')
load('/home/saintsbury/data/Dropbox/Projects/CVPR2022/Datasets/KIT_dataset/KIT_new2.mat')

% sampling parameters
point = point_no; % _super, _no; 104

rng(10);
color = rand(100, 3);
% color(4 : 8, :) = [];

%% prepare for data
rescale = 1;
downsample = 1;        
downsample_grid_size = 0.3; % 0.3
point_raw = cell(1, size(point, 2));
num_points = zeros(1, size(point, 2));
num_points_raw = zeros(1, size(point, 2));
for i = 1 : size(point, 2)
    if rescale == 1
        point{i} = point{i} - mean(point{i}, 2);
        max_length = max(max(point{1, i}));
        scale = max_length / 10;
        point{1, i} = point{1, i} / scale;
        point_raw{1, i} = point{1, i};
        num_points_raw(i) = size(point_raw{1, i}, 2);
    end
    if downsample == 1
        pc = pointCloud(point{1, i}');
        pc = pcdownsample(pc, 'gridAverage', downsample_grid_size);
        point{1, i} = pc.Location';
        num_points(i) = size(point{1, i}, 2);
        disp(['Case ', num2str(i), ' number of downsampled points = ', num2str(num_points(i))])
    else
        num_points(i) = size(point{1, i}, 2);
    end
end
disp('-------------------------------------------------------------------')
num_pc = size(point, 2);

%% select point cloud
% seal ---------------------------
% view_vector = [160 60];
% roll = 1600;
% turtle -------------------------
view_vector = [180 0];
roll = 0;
% cat1 ---------------------------
% view_vector = [120 0];
% roll = 90;
% cat2 dwarfs dog-----------------
% view_vector = [110 80];
% roll = 108;
% fish ---------------------------
% view_vector = [130 0];
% roll = 90;
% horse --------------------------
% view_vector = [90 90];
% roll = 90;
% bottles ------------------------
% view_vector = [30 -130];
% roll = -37;

close(figure(1))
figure(1)
j = 28; %17 28

showPoints(point{1, j}, 'MarkerSize', 10, 'ViewAxis', view_vector, 'CamRoll', roll)

%% fitting
para.w = 0.85;% 0.3 bottles \ 0.9 animals 0.8 seal 0.85
para.iterEM = 20; % 20
para.toleranceEM = 1e-3; % 1e-3
para.relative_toleranceEM = 0.2; % 0.2
para.iterLSQ_max = 2; % 2
para.sigma2 = 0.3; % 0.8 0.6 \ 0.3 0.2
para.max_switch = 2; % 2
para.adaptive_upper = 1; % 1
para.rescale = 0; % 1 \ 0

para.debug_mode = 0;
para.arclength = 0.2; %0.2

para.maximum_layer = 3; % 2 \ 3
para.minDistance = 1; % 1
para.minPoints = 60; % 60

[x, point_seg, point_outlier] = Hierarchical_SuperquadricsGaussian(point{j}, para);

%% fitting radial
% x = cell(1, 1);
% x{1} = cell(1, 1);
% x{1}{1} = superquadricsFitting(point{j}, 'Radial');

%% show segmentation
color_set = 1;
close(figure(2))
figure(2)
for i = 1 : size(point_seg, 2)
    for k = 1 : size(point_seg{i}, 2)
        hold on
        showPoints(point_seg{i}{k}, 'Color', color(color_set, :), 'MarkerSize', 80)       
        color_set = color_set + 1;
    end
end
hold off
view(view_vector)
camroll(roll)
%% show superquadrics
show_original_point = 0;

close(figure(3))
figure(3)
hold on
if show_original_point == 1
    showPoints(point{1, j}, 'MarkerSize', 10)
end
color_set = 1;

for i = 1 : size(x, 2)
    for k = 1 : size(x{i}, 2)
        hold on
        showSuperquadrics(x{i}{k}, 'Color', color(color_set, :), 'FaceAlpha', 1, 'Arclength', 0.1, 'Light', 0); %Aplha 0.6
%         showSuperquadrics(x{i}{k}, 'Color', [0 0 1], 'FaceAlpha', 1, 'Arclength', 0.1, 'Light', 0) %Aplha 0.6
        color_set = color_set + 1;
    end
end
hold off

view(view_vector)
camroll(roll)
light
material dull

%% show superquadrics
show_original_point = 0;

close(figure(4))
figure(4)
hold on
if show_original_point == 1
    showPoints(point{1, j}, 'MarkerSize', 10)
end
color_set = 1;

for i = 1 : size(x, 2)
    for k = 1 : size(x{i}, 2)
        hold on
        showSuperquadrics(x{i}{k}, 'Color', [0 0 1], 'FaceAlpha', 1, 'Arclength', 0.1, 'Light', 0); %Aplha 0.6
        color_set = color_set + 1;
    end
end
hold off

view(view_vector)
camroll(roll)
light
material dull