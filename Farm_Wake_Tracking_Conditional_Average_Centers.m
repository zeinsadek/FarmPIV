%% Quick proof-of-concept at phase-averaging based on wake position


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
tracking_folder = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/tracking';
experiment  = 'SingleFarm';
plane = 'Plane_1';


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD PIV RECORDINGS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find how many recordings are available in directory
files = dir(fullfile(data_folder, experiment));
names = {files.name};
names = erase(names, '.mat');
names = names(~contains(names, {'._', '..', '.'}));
recordings = sum(contains(names, plane));

% Load PIV data
for r = 1:recordings
    fprintf("Importing Recording %1.0f\n\n", r);
    name = strcat(plane, "_Recording_", num2str(r), "_DATA.mat");
    tmp = load(fullfile(data_folder, experiment, name));
    data(r) = tmp.output;
end

clear files names tmp name r recordings
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD WAKE CENTER TRACKING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find how many recordings are available in directory
files = dir(fullfile(tracking_folder, experiment));
names = {files.name};
names = erase(names, '.mat');
names = names(~contains(names, {'._', '..', '.'}));
recordings = sum(contains(names, plane));

% Load wake center data
for r = 1:recordings
    fprintf("Importing Recording %1.0f\n\n", r);
    name = strcat(plane, "_Recording_", num2str(r), "_WAKE_TRACKING.mat");
    tmp = load(fullfile(tracking_folder, experiment, name));
    tracking(r) = tmp.output;
end

clear files names tmp name r
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD ENSEMBLE AVERAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ensemble = load(fullfile(project_folder, "combined", experiment, strcat(plane, '_COMBINED.mat')));
ensemble = ensemble.output;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROP DATA DOWN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coordinates
D = 80;
Sy = 3;
X = (data(1).X + 200);
Y = (data(1).Y - 940);

% Crop
left_crop   = 1;
right_crop  = 4;
top_crop    = -12;
bottom_crop = 4.5;

[~, left_crop_index]   = min(abs(X(1,:) / D - left_crop));
[~, right_crop_index]  = min(abs(X(1,:) / D - right_crop));
[~, top_crop_index]    = min(abs(Y(:,1) / D - top_crop));
[~, bottom_crop_index] = min(abs(Y(:,1) / D - bottom_crop));

% Crop instantaneous
for r = 1:3
    data(r).U = data(r).U(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index, :);
    data(r).V = data(r).V(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index, :);
    clear r
end

X = X(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index);
Y = Y(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index);
x = X(1,:);
y = Y(:,1);

% Crop ensemble
components = {'u', 'v', 'uu', 'vv', 'uv'};
for c = 1:length(components)
    component = components{c};
    tmp = ensemble.(component);
    tmp = tmp(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index);
    ensemble.(component) = tmp;
    clear c tmp
end


figure()
contourf(X / D, Y / D, data(1).U(:,:,1), 100, 'linestyle', 'none')
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

clear left_crop right_crop top_crop left_crop_index right_crop_index top_crop_index 
clear bottom_crop bottom_crop_index components

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WAKE TRAJECTORY BINNING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Double check FPS
num_images = data(1).D;
PIV_FPS = 1385;
dt = 1/PIV_FPS;
time = linspace(0, dt * num_images, num_images);

% What angles to look for
target_angles = -15:15:15;

% What allowable tolerance for binning +/- [deg]
tolerance = 5;

% Colors
colors = parula(length(target_angles));


% Loop through conditioning different turbines
clc; close all
for reference_turbine = 1:5


    fprintf('Conditioning Turbine %1.0f\n\n', reference_turbine)
    for recording = 1:3
    
        % Which recording
        fprintf('Recording %1.0f:\n', recording)
    
        % Load signals for reference turbine
        reference_turbine_centers = tracking(recording).raw(reference_turbine).center_minimum;
        
        % Fit a line to each trajectory
        slopes = nan(size(reference_turbine_centers,1), 1);   
        for i = 1:num_images
            coefficients = polyfit(x, reference_turbine_centers(i,:), 1);
            slopes(i) = coefficients(1);
            clear i
        end
    
        % Convert slopes into angles [degrees]
        wake_angles = atan(slopes) * (180 / pi);

        % Bin each target angle
        for angle_enumerate = 1:length(target_angles)
        
            % Index target slope
            target_angle = target_angles(angle_enumerate);
            
            % Get upper and lower bounds
            lower_bound = target_angle - tolerance;
            upper_bound = target_angle + tolerance;
            fprintf('Lower bound: %2.1f, Upper bound: %2.1f\n', lower_bound, upper_bound)
        
            % Mask to find frames that fit the bounds
            mask = wake_angles > lower_bound & wake_angles < upper_bound;
            detected_images = sum(mask);
            fprintf('%3.0f Images fit\n\n', detected_images)

            % Save masks
            masks(recording, angle_enumerate).phase = mask;
    
            clear angle_enumerate
        end
    
        fprintf('\n')
        clear recording
    end

    
    % Print number of total images used in each phase
    for angle_enumerate = 1:length(target_angles)
        target_angle = target_angles(angle_enumerate);
    
        total_num_images = sum([masks(1, angle_enumerate).phase, masks(2, angle_enumerate).phase, masks(3, angle_enumerate).phase], 'all');
    
        fprintf('%2.1f wake angle uses %4.0f images\n', target_angle, total_num_images)
        clear angle_enumerate
    end
    
    clear coefficients detected_images dt h lower_bound upper_bound mask reference_turbine_centers
    clear wake_angles slopes target_angle total_num_images


    %%% Compute phase average
    % Loop through target angles
    for angle_enumerate = 1:length(target_angles)
        % Combine masked images across all recordings
        concatenated(angle_enumerate).U = cat(3, data(1).U(:,:,masks(1, angle_enumerate).phase), data(2).U(:,:,masks(2, angle_enumerate).phase), data(3).U(:,:,masks(3, angle_enumerate).phase));
        concatenated(angle_enumerate).V = cat(3, data(1).V(:,:,masks(1, angle_enumerate).phase), data(2).V(:,:,masks(2, angle_enumerate).phase), data(3).V(:,:,masks(3, angle_enumerate).phase));
        
        % Compute phase average
        phase_average(angle_enumerate) = data2means(concatenated(angle_enumerate));
    
        clear angle_enumerate
    end



    %%% Center line detection for phase average
    for angle_enumerate = 1:length(target_angles)
        phase_averaged_u = phase_average(angle_enumerate).u;
    
        % Find center line
        for turbine = 1:5
            % Crop just to turbine
            turbine_top_bound = -10.5 + Sy * (turbine - 1);
            turbine_bottom_bound = turbine_top_bound + Sy;
    
            fprintf('Turbine %1.0f: [%2.1f, %2.1f]\n', turbine, turbine_top_bound, turbine_bottom_bound)
    
            % Find indicies to crop to
            [~, turbine_top_bound_index] = min(abs(Y(:,1) / D - turbine_top_bound));
            [~, turbine_bottom_bound_index] = min(abs(Y(:,1) / D - turbine_bottom_bound));
            cropped_snapshot = phase_averaged_u(turbine_top_bound_index:turbine_bottom_bound_index, :);
            cropped_X = X(turbine_top_bound_index:turbine_bottom_bound_index, :);
            cropped_Y = Y(turbine_top_bound_index:turbine_bottom_bound_index, :);
            cropped_y = cropped_Y(:,1);
    
            %%% Part 1: Track minimum velocity
            [tst, min_velocity_indicies] = min(cropped_snapshot, [], 1);
    
            % Save to plot later
            % conditional_average_centers(turbine).center = (cropped_y(min_velocity_indicies) / D) - (-12 + turbine * 3);
            conditional_average_centers(angle_enumerate).angle(reference_turbine).constrained(turbine).center = (cropped_y(min_velocity_indicies) / D) - (-12 + turbine * 3);
        end
        clear angle_enumerate
    end
end
clear target_angle


%% 

lw = 3;
gaussian_window = 10;
turbine_colors = {'#F42C04', '#FFC857', '#018E42', '#3083DC', '#151E3F'};

clc; close all
figure('color', 'white')
tiledlayout(5,3)

c = 1;
for reference_turbine = 1:5
    for angle_enumerate = 1:length(target_angles)
        h(c) = nexttile;
        hold on
        for turbine = 1:5
            plot(x / D, smoothdata(conditional_average_centers(angle_enumerate).angle(reference_turbine).constrained(turbine).center, 'gaussian', 10), ...
                 'linewidth', lw, 'color', turbine_colors{turbine}, 'displayname', sprintf('Turbine %1.0f', turbine))
            scatter(x / D, smoothdata(conditional_average_centers(angle_enumerate).angle(reference_turbine).constrained(turbine).center, 'gaussian', 10), ...
                    20, 'filled', 'markerfacecolor', turbine_colors{turbine}, 'HandleVisibility', 'off')
        end
        hold off
        title(sprintf('Turbine %1.0f', reference_turbine))
        axis equal
        xlim([1, 4])
        ylim([-1.5, 1.5])
        set(gca, 'YDir', 'reverse')
        c = c + 1;
    end

    % legend('location', 'northeastoutside')
end
linkaxes(h, 'xy')





















%% Functions

function output = data2means(inst_struct)

    % Extract Instantaneous Velocities from Struct.
    inst_u  = inst_struct.U;
    inst_v  = inst_struct.V;


    % Calculate Velocity Means
    output.u = mean(inst_u, 3, 'omitnan');
    output.v = mean(inst_v, 3, 'omitnan');

    % Create Reynolds Stress Objects
    uu_p = zeros(size(inst_u));
    vv_p = zeros(size(inst_u));
    uv_p = zeros(size(inst_u));

    % Num frames
    D = size(inst_u, 3);

    % Loop Through Each Frame in Struct.
    fprintf('\n<data2means> PROGRESS: ');
    for frame_number = 1:D
        
        % Print Progress.
        progressbarText(frame_number/D);
        
        % Instantaneous Fluctuations.
        u_pi = inst_u(:, :, frame_number) - output.u;
        v_pi = inst_v(:, :, frame_number) - output.v;

        % Instantaneous Stresses.
        uu_pi = u_pi.*u_pi;
        vv_pi = v_pi.*v_pi;
        uv_pi = u_pi.*v_pi;

        % Array of Mean Stresses.
        uu_p(:, :, frame_number) = uu_pi;
        vv_p(:, :, frame_number) = vv_pi;
        uv_p(:, :, frame_number) = uv_pi;

    end
    
    % Mean Stresses.
    output.uu = mean(uu_p, 3, 'omitnan');
    output.vv = mean(vv_p, 3, 'omitnan');
    output.uv = mean(uv_p, 3, 'omitnan');
    output.D = D;

end
