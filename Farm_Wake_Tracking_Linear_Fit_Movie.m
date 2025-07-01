% Testing to see if a linear fit to the wake center is an appropriate
% simplification

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
project_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm';
data_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/data';
experiment  = 'SingleFarm';
recording   = 'Plane_1_Recording_1';

% Load PIV
data_path = fullfile(data_folder, experiment, strcat(recording, '_DATA.mat'));
data      = load(data_path);
data      = data.output;

% Load wake tracking
tracking = load(fullfile(project_folder, 'tracking', experiment, strcat(recording, '_WAKE_TRACKING.mat')));
tracking = tracking.output;

% Test method by looking at single turbine
u_inf  = 7.5;
Sy     = 3;

clear data_folder data_path

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

figure('color', 'white')
contourf(X / D, Y / D, inpaint_nans(U(:,:,1)), 100, 'linestyle', 'none')
axis equal
yticks(-12:3:3)
colorbar()
set(gca, 'YDir', 'reverse')
yline(4.5)
yline(1.5)
yline(-1.5)
yline(-4.5)
yline(-7.5)
yline(-10.5)

clear left_crop right_crop top_crop left_crop_index right_crop_index top_crop_index bottom_crop bottom_crop_index



%% Plot single instantaneous frame, wake centers, and linear fits

x = X(1,:);
frame = 1000;

figure('color', 'white')
hold on

% Plot velocity
contourf(X / D, Y / D, U(:,:,frame), 100, 'linestyle', 'none')

% Loop through turbine centers
for turbine = 1:5
    % Center shift
    center_shift = -9 + Sy * (turbine - 1);

    % Raw signal
    plot(X(1,:) / D, tracking.raw(turbine).center_minimum(frame, :) / D, 'color', 'black', 'linewidth', 2);
    % Filtered signal
    plot(x / D, tracking.raw(turbine).center_minimum(frame, :) / D, 'color', 'red', 'linewidth', 2);

    % Fit a line to the raw signal
    coefficients = polyfit(x / D, tracking.raw(turbine).center_minimum(frame, :) / D, 1);
    % Plot
    P = plot(x / D, polyval(coefficients, x / D), 'color', 'blue', 'linewidth', 2);
    P.Color(4) = 0.5;
    % Plot only with slope
    slope = coefficients(1);
    plot(x / D, slope * (x / D) - (slope) + center_shift, 'color', 'blue', 'linewidth', 2, 'LineStyle', '--')
end

% Axes size
axis equal
yticks(-12:3:3)
colorbar()
set(gca, 'YDir', 'reverse')



%% Loop and save as movie


% Video settings
FPS = 10;
skip = 5;
v = VideoWriter("movies/Wake_Tracking_Test_Linear_Fit_Filtered_Fit.mp4",'MPEG-4');
v.FrameRate = FPS;
open(v)


% Frame counter
c = 1;
num_images = 1000;

clc;
% Loop through frames
for frame = 1:skip:num_images

    % Progress bar
    progressbarText(c / (num_images /skip))

    % Start figure
    ax = figure('color', 'white', 'visible', 'off');
    hold on
    
    % Plot velocity
    contourf(X / D, Y / D, U(:,:,frame), 100, 'linestyle', 'none')
    
    % Loop through turbine centers
    for turbine = 1:5
        % Center shift
        center_shift = -9 + Sy * (turbine - 1);
    
        % Raw signal
        plot(X(1,:) / D, tracking.raw(turbine).center_minimum(frame, :) / D, 'color', 'black', 'linewidth', 2);
        % Filtered signal
        plot(x / D, tracking.filtered(turbine).center_minimum(frame, :) / D, 'color', 'blue', 'linewidth', 2);
    
        % Fit a line to the raw signal
        coefficients = polyfit(x / D, tracking.filtered(turbine).center_minimum(frame, :) / D, 1);
        % Plot
        P = plot(x / D, polyval(coefficients, x / D), 'color', 'red', 'linewidth', 2);
        P.Color(4) = 0.5;
        % Plot only with slope
        slope = coefficients(1);
        plot(x / D, slope * (x / D) - (slope) + center_shift, 'color', 'red', 'linewidth', 2, 'LineStyle', '--')
    end
    
    % Axes size
    axis equal
    yticks(-12:3:3)
    colorbar()
    set(gca, 'YDir', 'reverse')
    clim([0.2, 1])

    % Save to video
    videoFrame = getframe(ax);
    close all
    writeVideo(v, videoFrame);

    % Increment counter
    c = c + 1;
end

% Stop video
close(v)
clc;








