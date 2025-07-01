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
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Data paths
data_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/data';
experiment  = 'SingleFarm';
% experiment = 'Farm2Farm_40D_Gap';
recording   = 'Plane_1_Recording_1';

data_path = fullfile(data_folder, experiment, strcat(recording, '_DATA.mat'));
data      = load(data_path);
data      = data.output;

% Test method by looking at single turbine
u_inf  = 8;
Sy     = 3;

clear data_folder recording data_path

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROP DATA DOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates
D = 80;
X = (data.X + 200);
Y = (data.Y - 940);

% Instantaneous snapshots
U = double(data.U) / u_inf;

% Crop
left_crop   = 1;
right_crop  = 4;
top_crop    = -12;
bottom_crop = 4.5;

[~, left_crop_index]   = min(abs(X(1,:) / D - left_crop));
[~, right_crop_index]  = min(abs(X(1,:) / D - right_crop));
[~, top_crop_index]    = min(abs(Y(:,1) / D - top_crop));
[~, bottom_crop_index] = min(abs(Y(:,1) / D - bottom_crop));

U = U(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index, :);
X = X(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index);
Y = Y(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index);

figure()
contourf(X/D, Y/D, U(:,:,1), 100, 'linestyle', 'none')
axis equal
yticks(-12:3:3)
colorbar()
set(gca, 'YDir', 'reverse')
yline(1.5)
yline(-1.5)
yline(-4.5)
yline(-7.5)
yline(-10.5)

clear left_crop right_crop top_crop left_crop_index right_crop_index top_crop_index

%% See if we can get a line for each wake from a single snapshot


snapshot = U(:,:,1);

% Fill in any nans
snapshot = inpaint_nans(snapshot);

% Crop to a single turbine for now: center one
turbine_top_bound = -10.4;
turbine_bottom_bound = -7.5;

[~, turbine_top_bound_index] = min(abs(Y(:,1) / D - turbine_top_bound));
[~, turbine_bottom_bound_index] = min(abs(Y(:,1) / D - turbine_bottom_bound));

snapshot = snapshot(turbine_top_bound_index:turbine_bottom_bound_index, :);
cropped_X = X(turbine_top_bound_index:turbine_bottom_bound_index, :);
cropped_Y = Y(turbine_top_bound_index:turbine_bottom_bound_index, :);

%% Try blurring it a little
snapshot_blurred = imgaussfilt(snapshot, 1);

% Plot to compare
figure()
tiledlayout(1,2)
nexttile
contourf(cropped_X / D, cropped_Y / D, snapshot, 100, 'linestyle', 'none')
axis equal
set(gca, 'YDir', 'reverse')
title('Original')

nexttile
contourf(cropped_X / D, cropped_Y / D, snapshot_blurred, 100, 'linestyle', 'none')
axis equal
set(gca, 'YDir', 'reverse')
title('Blurred')


%% Method 1: Track minimum velocity 

% Using blurred snapshot to determine centerline
[~, min_velocity_indicies] = min(snapshot_blurred, [], 1);

figure()
tiledlayout(1,2)
nexttile
hold on
contourf(cropped_X / D, cropped_Y / D, snapshot, 100, 'linestyle', 'none')
plot(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 'color', 'red', 'linewidth', 2)
scatter(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 30, 'filled', 'MarkerFaceColor', 'red')
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('Original')

nexttile
hold on
contourf(cropped_X / D, cropped_Y / D, snapshot_blurred, 100, 'linestyle', 'none')
plot(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 'color', 'red', 'linewidth', 2)
scatter(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 30, 'filled', 'MarkerFaceColor', 'red')
hold off
axis equal
set(gca, 'YDir', 'reverse')
title('Blurred')


%% Method 2: Gaussian fitting

gaussian_centers = nan(1, size(snapshot, 2));
gaussian_widths  = nan(1, size(snapshot, 2));
% fit_success = false(1, nx);

% Gaussian equation and initial guesses
gaussEqn = @(b, y) b(1) * exp(-((y - b(2)).^2) / (2 * b(3)^2));
y = cropped_Y(:,1);
S0 = 40;

for s = 1:size(snapshot, 2)

    % Compute velocity deficit
    deficit = 1 - snapshot_blurred(:, s);
    
    % Initial guesses for amplitude and center based on data and minimum
    % velocity fitting
    A0 = max(deficit, [], 'all');
    y0 = y(min_velocity_indicies(s));
    
    % Fit using nonlinear least squares
    try
        opts = optimoptions('lsqcurvefit', ...
        'FunctionTolerance', 1e-10, ...
        'MaxIterations', 1000, ...
        'Display', 'off');
        beta = lsqcurvefit(gaussEqn, [A0, y0, S0], y, deficit, [], [], opts);
        gaussian_centers(s) = beta(2);
        gaussian_widths(s)  = beta(3);
    catch
        % If fitting fails, leave as NaN
        fprintf('Failure!\n')
    end
end


figure()
hold on
% Instantaneous data
contourf(cropped_X / D, cropped_Y / D, snapshot, 100, 'linestyle', 'none')

% Minimum velocity-tracled center
plot(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 'color', 'red', 'linewidth', 2)
scatter(cropped_X(1,:) / D, cropped_Y(min_velocity_indicies, :) / D, 30, 'filled', 'MarkerFaceColor', 'red')

% Gaussian-tracked center
plot(cropped_X(1,:) / D, gaussian_centers / D, 'color', 'blue', 'linewidth', 2)
scatter(cropped_X(1,:) / D, gaussian_centers / D, 30, 'filled', 'MarkerFaceColor', 'blue')

% Gaussian width
plot(cropped_X(1,:) / D, (gaussian_centers / D) + (gaussian_widths / D), 'color', 'green', 'linewidth', 2)
plot(cropped_X(1,:) / D, (gaussian_centers / D) - (gaussian_widths / D), 'color', 'green', 'linewidth', 2)

hold off
axis equal
set(gca, 'YDir', 'reverse')
title('Original')















