%%% Wake Tracking
% Zein Sadek
% PSU + Oldenburg

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions/Inpaint_nans/Inpaint_nans');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/colormaps');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Data paths
data_folder = '/Volumes/Zein_PIV_4/Oldenburg_Convergence/data';
recording   = 'Plane_1_Recording_1';

data_path   = fullfile(data_folder, strcat(recording, '_DATA.mat'));
data        = load(data_path);
data        = data.output;

% Test method by looking at single turbine
u_inf  = 8;     % m/s
offset = -208;  % mm
D      = 80;    % mm
Sy     = 3;     % D

% Coordinates
X = data.X.';
Y = data.Y.';
Y = Y - offset;

% Instantaneous snapshots
U = data.U;
V = data.V;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAKE CENTER TRACKING (TESTING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for j = 1:10

    % Index single image
    snapshot = U(:,:,j);
    snapshot = snapshot / u_inf;
    
    % Crop to single turbine
    snapshot(Y > (Sy * D) / 2 | Y < (-Sy * D) / 2) = nan;
    
    % Finding location of minimum velocity within profiles
    y           = Y(1,:);
    profile     = snapshot(1,:);
    tmp         = size(snapshot);
    wake_center = nan(1,tmp(1));
    
    % Iterate through profiles
    for i = 1:tmp(1)
        profile = snapshot(i,:);
        % Find minimum/peak of profile
        [M,I] = min(profile, [], 'all');
        wake_center(i) = y(I);
    end
    
    % Lowpass filter to smooth
    wake_center_lowpass = lowpass(wake_center, 0.05);
    
    % Plot
    figure()
    colormap(parula)
    hold on
    
    % Snapshot
    contourf(X,Y,snapshot, 500, linestyle = 'none')
    
    % Unfiltered Wake Center
    plot(unique(X), wake_center, 'r')
    scatter(unique(X), wake_center, 10, 'filled', 'MarkerFaceColor', 'r')
    
    % Lowpass Filtered Wake Center
    plot(unique(X), wake_center_lowpass, 'w')
    scatter(unique(X), wake_center_lowpass, 10, 'filled', 'MarkerFaceColor', 'w')
    
    hold off
    axis equal
    colorbar
    xlim([min(X, [], 'all'),max(X, [], 'all')])
    ylim([-Sy * D / 2, Sy * D / 2])
end

% %% Wake Edge
% 
% copy = imgaussfilt(snapshot, 2);
% threshold = 0.85;
% copy(copy >= threshold) = nan;
% % copy(copy < threshold) = 1;
% 
% figure()
% contourf(X,Y,copy, 500, linestyle = 'none')
% axis equal
% colorbar
% xlim([min(X, [], 'all'),max(X, [], 'all')])
% ylim([-Sy * D / 2, Sy * D / 2])











%% Gaussian Fitting
% clc;
% y = Y(1,:);
% % Fitting gaussian curve to profile
% profile = snapshot(1,:);
% gaussian = @(x, ypos) x(1) * (exp(-(1/2) * ((ypos - x(2)) / x(3)).^2)) + x(4);
% x0 = [1.1, 0, 35, 0];
% 
% nan_mask = ~isnan(profile);
% options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt', 'MaxIterations', 1E8);
% tst = lsqcurvefit(gaussian, x0, double(profile(nan_mask)), double(y(nan_mask)), [-2, -50, 0, 0], [0, 50, 100, 2], options)
% 
% figure()
% hold on
% plot(y, -profile + 1)
% plot(y(nan_mask), gaussian(x0, y(nan_mask)))
% hold off












