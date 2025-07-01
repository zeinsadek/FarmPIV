%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps');
fprintf("All Paths Imported...\n\n");
% addpath('G:/Other computers/Zein MacBook Pro/Farm/Farm_Functions/');
% addpath('C:/Users/Zein/Documents/MATLAB/readimx-v2.1.8-win64/');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data paths
experiment = 'Farm2Farm_40D_Gap';
recording_name = 'Plane_9_Recording_3';

project_path   = fullfile('/Volumes/Frm2FrmProc/', experiment, 'Oldenburg_Farm2Farm_40D_Gap_Block3');
processing     = 'TR_PIV_MPd(2x32x32_50%ov_ImgCorr)_GPU';
inpt_name      = recording_name;


% Image paths
piv_path = fullfile(project_path, recording_name, 'ImageCorrection', processing);
% piv_path = fullfile(project_path, recording_name, processing);

% Save paths
results_path = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm';
mtlb_file    = fullfile(results_path, 'data'   , experiment, strcat(inpt_name, '_DATA.mat'));
mean_file    = fullfile(results_path, 'means'  , experiment, strcat(inpt_name, '_MEANS.mat'));
figure_file  = fullfile(results_path, 'figures', experiment);

if ~exist(fullfile(results_path, 'data', experiment), 'dir')
    mkdir(fullfile(results_path, 'data', experiment))
end

if ~exist(fullfile(results_path, 'means', experiment), 'dir')
    mkdir(fullfile(results_path, 'means', experiment))
end

if ~exist(fullfile(results_path, 'figures', experiment), 'dir')
    mkdir(fullfile(results_path, 'figures', experiment))
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DAVIS TO MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mtlb_file, 'file')
    fprintf('* Loading DATA from File\n')
    data = load(mtlb_file);
    data = data.output;
else
    tic
    data = vector2matlab(piv_path, mtlb_file);
    toc
end

% data = vector2matlab(piv_path, mtlb_file);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB DATA TO ENSEMBLE MEANS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(mean_file, 'file')
    fprintf('* Loading MEANS from File\n')
    means = load(mean_file); 
    means = means.output;
else
    tic
    means = data2means(mean_file, data);
    toc
end

% means = data2means(mean_file, data);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;

% Coordinates
D = 80;
X = (means.X + 200) / D;
Y = (means.Y - 940) / D;

% Means
U = (means.u);
V = (means.v);

% Stresses
uu = means.uu;
vv = means.vv;
uv = means.uv;

% Means Plots

xLower = -14;
xUpper = 4.5;

close all
ax = figure();
t  = tiledlayout(1,2);
sgtitle(inpt_name, 'interpreter', 'none')

ax1 = nexttile();
hold on
colormap(ax1, jet)
contourf(X, Y, U, 500, 'linestyle', 'none')
axis equal
ylim([xLower,xUpper])
xlim([1,4])
yticks(-12:3:xUpper)
clim([1, 8])
colorbar()
title('u')
hold off

ax2 = nexttile();
hold on
colormap(ax2, coolwarm);
contourf(X, Y, V, 500, 'linestyle', 'none')
axis equal
ylim([xLower,xUpper])
xlim([1,4])
yticks(-12:3:xUpper)
clim([-0.4, 0.4])
colorbar()
title('v')
hold off

set([ax1, ax2], 'YDir','reverse')

% exportgraphics(t, fullfile(figure_file, strcat(recording_name, '_U_V.png')), 'resolution', 300)



%% Stresses Plots

ax = figure();
t  = tiledlayout(1,3);
sgtitle(inpt_name, 'interpreter', 'none')

% Normal Stresses
nexttile()
colormap jet
contourf(X, Y, uu, 500, 'linestyle', 'none')
axis equal
ylim([xLower,xUpper])
xlim([1,4])
yticks(-3:3:xUpper)
colorbar()
title('uu')

nexttile()
colormap jet
contourf(X, Y, vv, 500, 'linestyle', 'none')
axis equal
ylim([xLower,xUpper])
xlim([1,4])
yticks(-3:3:xUpper)
colorbar()
title('vv')

% Shear Stresses
nexttile()
colormap jet
contourf(X, Y, uv, 500, 'linestyle', 'none')
axis equal
ylim([xLower,xUpper])
xlim([1,4])
yticks(-3:3:xUpper)
colorbar()
title('uv')

%% Convergence in Means

% clc;
% N      = data.D;
% step   = 50;
% groups = 1:step:N;
% avgs_wake  = zeros(size(groups));
% avgs_aisle = zeros(size(groups));
% avgs_edge  = zeros(size(groups));
% avgs_all   = zeros(size(groups));
% 
% c = 1;
% for i = step:step:N
%     disp(i)
%     % Wake
%     avgs_wake(c) = mean(data.U(r, c_center_wake, 1:i), 'omitnan');
% 
%     % Aisle
%     avgs_aisle(c) = mean(data.U(r, c_aisle, 1:i), 'omitnan');
% 
%     % Edge
%     avgs_edge(c) = mean(data.U(r, c_edge, 1:i), 'omitnan');
% 
%     % All
%     avgs_all(c) = mean(data.U(:, :, 1:i), 'all', 'omitnan');
%     c = c + 1;
% end
% 
% %% Plot Convergence
% 
% ls = 2;
% figure();
% hold on
% plot(groups, avgs_wake, 'displayname', 'Wake', 'linewidth', ls)
% plot(groups, avgs_aisle, 'displayname', 'Aisle', 'linewidth', ls)
% plot(groups, avgs_edge, 'displayname', 'Edge', 'linewidth', ls)
% plot(groups, avgs_all, 'displayname', 'All', 'linewidth', ls)
% hold off
% 
% legend()
% xlim([step, means.D - 2 * step]);
% ylim([0, 10]);
% xlabel('Num Images');
% ylabel('Mean [m/s]');
% title('Convergence')
















