%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/colormaps');
fprintf("All Paths Imported...\n\n");
% addpath('G:/Other computers/Zein MacBook Pro/Farm/Farm_Functions/');
% addpath('C:/Users/Zein/Documents/MATLAB/readimx-v2.1.8-win64/');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
project_path   = '/Volumes/Farm/PIV/PSU_Farm_Single_Farm_Date=20240811';
recording_name = 'Plane_1_Recording_3';
processing     = 'TR_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU';
inpt_name      = recording_name;


% Image paths
piv_path = fullfile(project_path, recording_name, processing);

% Save paths
results_path = '/Volumes/ATHENA/TurbulenceClass';
mtlb_file    = fullfile(results_path, 'data'   , strcat(inpt_name, '_DATA.mat'));
mean_file    = fullfile(results_path, 'means'  , strcat(inpt_name, '_MEANS.mat'));
figure_file  = fullfile(results_path, 'figures');

% Make specific folder for figures of an experiment
% if ~exist(figure_file, 'dir')
%     mkdir(figure_file)
% end

if ~exist(fullfile(results_path, 'data'), 'dir')
    mkdir(fullfile(results_path, 'data'))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAVIS TO MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mtlb_file, 'file')
    fprintf('* Loading DATA from File\n')
    data = load(mtlb_file);
    data = data.output;
else
    data = vector2matlab(piv_path, mtlb_file);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB DATA TO ENSEMBLE/PHASE MEANS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mean_file, 'file')
     fprintf('* Loading MEANS from File\n')
     means = load(mean_file); 
     means = means.output;
else
     means = data2means(mean_file, data);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

% Coordinates
X = means.X.';
Y = means.Y.';

% Zero to center turbine
% Y = Y - 1125;

% Means
U = means.u;
V = means.v;

% Stresses
uu = means.uu;
vv = means.vv;
uv = means.uv;

%% Means Plots

r = 20;
c_center_wake = 403;
c_aisle = 307;
c_edge = 130;

ax = figure();
t  = tiledlayout(1,2);
sgtitle(inpt_name, 'interpreter', 'none')

ax1 = nexttile();
hold on
colormap(ax1, jet)
contourf(X, Y, U, 500, 'linestyle', 'none')
% scatter3(X(r,c_center_wake), Y(r,c_center_wake), 0, 50, 'k', 'filled')
% scatter3(X(r,c_aisle), Y(r,c_aisle), 0, 50, 'k', 'filled')
% scatter3(X(r,c_edge), Y(r,c_edge), 0, 50, 'k', 'filled')
axis equal
xlim([-120,120])
ylim([-1300,150])
caxis([-1, 8])
colorbar()
title('u')
hold off

ax2 = nexttile();
hold on
colormap(ax2, coolwarm);
contourf(X, Y, V, 500, 'linestyle', 'none')
axis equal
xlim([-120,120])
ylim([-1300,150])
caxis([-0.4, 0.4])
colorbar()
title('v')
hold off

exportgraphics(t, fullfile(figure_file, strcat(recording_name, '_U_V.png')), 'resolution', 300)



%% Stresses Plots

ax = figure();
t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

% Normal Stresses
nexttile()
colormap jet
contourf(X, Y, uu, 500, 'linestyle', 'none')
axis equal
xlim([-120,120])
ylim([-1300,150])
colorbar()
title('uu')

nexttile()
colormap jet
contourf(X, Y, vv, 500, 'linestyle', 'none')
axis equal
xlim([-120,120])
ylim([-1300,150])
colorbar()
title('vv')

% Shear Stresses
nexttile()
colormap jet
contourf(X, Y, uv, 500, 'linestyle', 'none')
axis equal
xlim([-120,120])
ylim([-1300,150])
colorbar()
title('uv')

%% Convergence in Means

clc;
N      = data.D;
step   = 50;
groups = 1:step:N;
avgs_wake  = zeros(size(groups));
avgs_aisle = zeros(size(groups));
avgs_edge  = zeros(size(groups));
avgs_all   = zeros(size(groups));

c = 1;
for i = step:step:N
    disp(i)
    % Wake
    avgs_wake(c) = mean(data.U(r, c_center_wake, 1:i), 'omitnan');
    
    % Aisle
    avgs_aisle(c) = mean(data.U(r, c_aisle, 1:i), 'omitnan');
    
    % Edge
    avgs_edge(c) = mean(data.U(r, c_edge, 1:i), 'omitnan');
    
    % All
    avgs_all(c) = mean(data.U(:, :, 1:i), 'all', 'omitnan');
    c = c + 1;
end

%% Plot Convergence

ls = 2;
figure();
hold on
plot(groups, avgs_wake, 'displayname', 'Wake', 'linewidth', ls)
plot(groups, avgs_aisle, 'displayname', 'Aisle', 'linewidth', ls)
plot(groups, avgs_edge, 'displayname', 'Edge', 'linewidth', ls)
plot(groups, avgs_all, 'displayname', 'All', 'linewidth', ls)
hold off

legend()
xlim([step, means.D - 2 * step]);
ylim([0, 10]);
xlabel('Num Images');
ylabel('Mean [m/s]');
title('Convergence')
















