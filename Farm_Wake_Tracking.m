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
project_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm';
data_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/data';
experiment  = 'SingleFarm';
recording   = 'Plane_2_Recording_3';

% Load PIV
data_path = fullfile(data_folder, experiment, strcat(recording, '_DATA.mat'));
data      = load(data_path);
data      = data.output;

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

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRACK CENTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Detection settings
 num_images = data.D;
 % num_images = 1000;

% smoothing_kernel = 10;
smoothing_kernel = 2;
frame_counter = 1;

% Initialize once before looping
for turbine = 1:5
    tracking(turbine).center_gaussian = nan(num_images, size(X, 2));
    tracking(turbine).center_minimum  = nan(num_images, size(X, 2));
    tracking(turbine).edge_gaussian   = nan(num_images, size(X, 2));
end

clc; close all;
% Loop through different frames
for frame = 1:num_images
    
    % Get a snapshot
    % fprintf("Frame %4.0f\n", frame)
    progressbarText(frame / num_images)
    snapshot = U(:,:,frame);
    
    % Fill in any nans
    snapshot = inpaint_nans(snapshot);
    
    % Loop through the different turbines
    for turbine = 1:5
    
        % Crop just to turbine
        % turbine_top_bound = 1.5 - (Sy * (turbine - 1));
        % turbine_bottom_bound = 4.5 - (Sy * (turbine - 1));
        turbine_top_bound = -10.5 + Sy * (turbine - 1);
        turbine_bottom_bound = turbine_top_bound + Sy;

        % fprintf('Turbine %1.0f: [%2.1f, %2.1f]\n', turbine, turbine_top_bound, turbine_bottom_bound)
        
        % Find indicies to crop to
        [~, turbine_top_bound_index] = min(abs(Y(:,1) / D - turbine_top_bound));
        [~, turbine_bottom_bound_index] = min(abs(Y(:,1) / D - turbine_bottom_bound));
        cropped_snapshot = snapshot(turbine_top_bound_index:turbine_bottom_bound_index, :);
        cropped_X = X(turbine_top_bound_index:turbine_bottom_bound_index, :);
        cropped_Y = Y(turbine_top_bound_index:turbine_bottom_bound_index, :);
          
        %%% Part 1: Track minimum velocity 
        % Using blurred snapshot to determine centerline

        % 2D Blur
        snapshot_blurred = imgaussfilt(cropped_snapshot, smoothing_kernel);

        %%% Part 1: Track minimum velocity
        [~, min_velocity_indicies] = min(snapshot_blurred, [], 1);
    
        %%% Part 2: Gaussian fitting
        gaussian_centers = nan(1, size(cropped_snapshot, 2));
        gaussian_widths  = nan(1, size(cropped_snapshot, 2));

        % Gaussian equation and initial guesses
        gaussEqn = @(b, y) b(1) * exp(-((y - b(2)).^2) / (2 * b(3)^2));
        y = cropped_Y(:,1);
        S0 = 40;
    
        % for slice = 1:size(cropped_snapshot, 2)
        % 
        %     % Compute velocity deficit
        %     deficit = 1 - snapshot_blurred(:, slice);
        % 
        %     % Initial guesses for amplitude and center based on data and minimum
        %     % velocity fitting
        %     A0 = max(deficit, [], 'all');
        %     y0 = y(min_velocity_indicies(slice));
        % 
        %     % Bounds for fitting
        %     lb = [0.2, min(y, [], 'all'), 0];
        %     ub = [1.0, max(y, [], 'all'), 200];
        % 
        %     % Fit using nonlinear least squares
        %     try
        %         opts = optimoptions('lsqcurvefit', ...
        %         'FunctionTolerance', 1e-10, ...
        %         'MaxIterations', 1000, ...
        %         'Display', 'off');
        %         beta = lsqcurvefit(gaussEqn, [A0, y0, S0], y, deficit, lb, ub, opts);
        %         gaussian_centers(slice) = beta(2);
        %         gaussian_widths(slice)  = beta(3);
        %     catch
        %         % If fitting fails, leave as NaN
        %         fprintf('Failure!\n')
        %     end
        % end


        % Save for all looped snapshots
        % tracking(turbine).center_gaussian(frame_counter, :) = gaussian_centers;
        tracking(turbine).center_minimum(frame_counter, :)  = y(min_velocity_indicies);
        % tracking(turbine).edge_gaussian(frame_counter, :)   = gaussian_widths;
  
    end

    % Increment frame counter
    frame_counter = frame_counter + 1;
end


clc;
clear frame_counter frame turbine slice A0 beta centers cropped_snapshot cropped_X cropped_Y deficit gaussian_centers 
clear gaussian_widths min_velocity_indicies opts s S0 snapshot snapshot_blurred turbine_bottom_bound_index turbine_bottom_bound
clear turbine_top_bound turbine_top_bound_index x_fill y0 y_fill lb ub gaussEqn


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FILTER DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filtered_tracking = tracking;
time_smoothing_kernel = 5;
gradient_threshold = 3.0;
x = X(1,:);

clc;
% Place nans where spikes occur
for turbine = 1:5
    for frame = 1:num_images
        % Load data
        center_gaussian = tracking(turbine).center_gaussian(frame,:);
        center_minimum = tracking(turbine).center_minimum(frame,:);

        % Detect and delete large spikes
        gaussian_spike_index = find(abs(gradient(center_gaussian, x)) > gradient_threshold, 1, 'first');
        minimum_spike_index = find(abs(gradient(center_minimum, x)) > gradient_threshold, 1, 'first');

        % Replace with nans
        center_gaussian(gaussian_spike_index:end) = nan;
        center_minimum(minimum_spike_index:end) = nan;

        % Save
        filtered_tracking(turbine).center_gaussian(frame, :) = center_gaussian;
        filtered_tracking(turbine).center_minimum(frame, :) = center_minimum;
    end
end


% Fill in placed nans
for turbine = 1:5

    % Fill in nans: spatio-temporal
    fprintf('Turbine %1.0f\n', turbine)

    % Cubic interpolation
    center_gaussian = fillmissing2(filtered_tracking(turbine).center_gaussian, 'cubic');
    center_minimum = fillmissing2(filtered_tracking(turbine).center_minimum, 'cubic');

    % Linear fallback
    center_gaussian = fillmissing2(center_gaussian, 'linear');
    center_minimum = fillmissing2(center_minimum, 'linear');

    % Nearest fallback
    center_gaussian = fillmissing2(center_gaussian, 'nearest');
    center_minimum = fillmissing2(center_minimum, 'nearest');

    % Smooth along columns
    % center_gaussian = smoothdata(center_gaussian, 1, 'gaussian', time_smoothing_kernel, 'omitnan');
    center_minimum = smoothdata(center_minimum, 1, 'gaussian', time_smoothing_kernel, 'omitnan');

    % Save back to array
    filtered_tracking(turbine).center_gaussian = center_gaussian;
    filtered_tracking(turbine).center_minimum = center_minimum;
end

clear turbine center_gaussian center_minimum frame gaussian_spike_index minimum_spike_index

%% Plot: time series

% Double check FPS
PIV_FPS = 1385;
dt = 1/PIV_FPS;
time = linspace(0, dt * num_images, num_images);

clc; close all; figure('color', 'white')
tiledlayout(5,1)

index = 50;
for turbine = 1:5
    h(turbine) = nexttile;

    % Center shift
    center_shift = -9 + Sy * (turbine - 1);

    % Raw: Minimum
    raw_centered_minimum = tracking(turbine).center_minimum(:,index)/D - center_shift;

    % Raw: Center data gaussian
    raw_centered_gaussian = tracking(turbine).center_gaussian(:,index)/D - center_shift;
        
    % Filtered: Center data gaussian
    filtered_centered_gaussian = filtered_tracking(turbine).center_gaussian(:,index)/D - center_shift;

    % Filtered: Center data minimum
    filtered_centered_minimum = filtered_tracking(turbine).center_minimum(:,index)/D - center_shift;

    hold on
    plot(time, raw_centered_minimum, 'color', 'blue')
    % plot(time, raw_centered_gaussian, 'color', 'red')
    plot(time, filtered_centered_gaussian, 'color', 'black')
    plot(time, filtered_centered_minimum, 'color', 'red')
    hold off
    title(sprintf('Turbine %1.0f', turbine))
    ylim([-Sy/2, Sy/2])
    xlim([0, max(time)])
    yline(0)
end


% linkaxes(h, 'xy')

clear dt filtered_centered_gaussian filtered_centered_minimum FPS h index PIV_FPS raw_centered_gaussian raw_centered_minimum time


%% Plot instantaneous centers

x = X(1,:);

colors = parula(num_images);
alpha = 1;
skip = 100;

% Plot all profiles for each turbuine
figure()
tiledlayout(5,1)

for turbine = 1:5
    h(turbine) = nexttile;
    hold on
    for frame = 1:skip:num_images

        % Center shift
        center_shift = -9 + Sy * (turbine - 1);

        % Raw: Center data
        raw_centered_data = tracking(turbine).center_minimum(frame,:)/D - center_shift;
        
        % Filtered: Center data gaussian
        % gaussian_filtered_centered_data = filtered_tracking(turbine).center_gaussian(frame,:)/D;% - center_shift;

        % Filtered: Center data minimum
        minimum_filtered_centered_data = filtered_tracking(turbine).center_minimum(frame,:)/D - center_shift;

        % Plot
        plot(x/D, raw_centered_data, 'color', 'red')
        % plot(x/D, gaussian_filtered_centered_data, 'color', 'black')
        plot(x/D, minimum_filtered_centered_data, 'color', 'blue')
    end
    hold off
    title(sprintf("Turbine %1.0f", turbine))
    yline(0)
    axis equal
    xlim([1,4])
    ylim([-Sy/2, Sy/2])
    ylabel('$z / D$', 'interpreter', 'latex')
    set(gca, 'YDir', 'reverse')
end

linkaxes(h, 'xy')

clear alpha colors frame gaussian_filtered_centered_data minimum_filtered_centered_data raw_centered_data skip turbine


%% Movie of filtered and raw data to visually compare

% Video settings
FPS = 10;
skip = 10;

% Save settings
movie_name = strcat(experiment, '_', recording, '_WakeCenter.mp4');
v = VideoWriter(fullfile('movies', movie_name),'MPEG-4');
v.FrameRate = FPS;
open(v)

sz = 5;

c = 1;
for frame = 1:skip:1000

    progressbarText(c / (num_images /skip))
    ax = figure('color', 'white', 'visible', 'off');
    hold on
    % contourf(X / D, Y / D, inpaint_nans(filloutliers(U(:,:,frame), nan)), 20, 'linestyle', 'none')
    contourf(X / D, Y / D, U(:,:,frame), 20, 'linestyle', 'none')

    for turbine = 1:5
        % Raw
        % plot(x / D, tracking(turbine).center_gaussian(frame, :) / D, 'color', 'red')
        % scatter(x / D, tracking(turbine).center_gaussian(frame, :) / D, sz, 'filled', 'MarkerFaceColor', 'red')

        % Raw: Minimum
        plot(x / D, tracking(turbine).center_minimum(frame, :) / D, 'color', 'red')
        scatter(x / D, tracking(turbine).center_minimum(frame, :) / D, sz, 'filled', 'MarkerFaceColor', 'red')

        % Filtered: Minimum
        plot(x / D, filtered_tracking(turbine).center_minimum(frame, :) / D, 'color', 'green')
        scatter(x / D, filtered_tracking(turbine).center_minimum(frame, :) / D, sz, 'filled', 'MarkerFaceColor', 'green')

        % Filtered: Gaussian
        % plot(x / D, filtered_tracking(turbine).center_gaussian(frame, :) / D, 'color', 'blue')
        % scatter(x / D, filtered_tracking(turbine).center_gaussian(frame, :) / D, sz, 'filled', 'MarkerFaceColor', 'blue')
    end

    hold off
    axis equal
    xlim([1,4]); ylim([-10.5, 4.5])
    colorbar(); colormap(bone); clim([0, 1])
    set(gca, 'YDir', 'reverse'); yticks(-12:3:3)
    yline(4.5)
    yline(1.5)
    yline(-1.5)
    yline(-4.5)
    yline(-7.5)
    yline(-10.5)

    % Save to video
    videoFrame = getframe(ax);
    close all
    writeVideo(v, videoFrame);

    c = c + 1;
end

close(v);

%% SAVE TO A MATFILE

output.raw = tracking;
output.filtered = filtered_tracking;
output.smoothing_kernel = smoothing_kernel;
output.x = X(1,:);

save_dir = fullfile(project_folder, 'tracking', experiment);
save_name = strcat(recording, '_WAKE_TRACKING.mat');
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

fprintf('Saving matfile...\n')
save(fullfile(save_dir, save_name), 'output')
clc; fprintf('Matfile saved!\n')



