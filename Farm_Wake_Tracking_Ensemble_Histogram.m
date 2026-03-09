%% Plot histogram of turbine wake centers to see if there are preferential directions


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
% FIT A LINE TO ALL WAKE CENTERS AND PLOT HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of images
num_images = length(tracking(1).raw(1).center_gaussian);

% Plot as a histogram
nbins = 9;
clc; close all
figure('color', 'white')
tiledlayout(5, 1, 'tilespacing', 'compact')


% Loop through turbines
for turbine = 1:5

    disp(turbine)
    h(turbine) = nexttile;
    ax = gca;

    % Fit a line to each trajectory
    slopes = nan(3 * num_images, 1);   
    
    for recording = 1:3
        for i = 1:num_images
            coefficients = polyfit(tracking(1).x, tracking(recording).raw(turbine).center_minimum(i,:), 1);
            slopes(i) = coefficients(1);
            clear i
        end
        clear recording
    end
    
    % Convert slopes into angles [degrees]
    wake_angles = atan(slopes) * (180 / pi);
    mean_wake_angle = mean(wake_angles, 'all', 'omitnan');

    if sign(mean_wake_angle) == 1
        color = 'green';
    elseif sign(mean_wake_angle) == -1
        color = 'red';
    else
        color = 'black';
    end

    hold on
    histogram(wake_angles, nbins, 'normalization', 'probability', 'orientation', 'horizontal'); %, 'FaceColor', color)
    % yline(mean_wake_angle)
    hold of
    title(sprintf('Turbine %1.0f', turbine))
    ylim([-30, 30])

    if turbine ~= 5
        ax.XTickLabel = [];
    end

    clear slopes wake_angles turbine clear ax
end

linkaxes(h, 'xy')



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIT A LINE TO ALL WAKE CENTERS AND PLOT POLAR HISTOGRAM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of images
num_images = length(tracking(1).raw(1).center_gaussian);

% Plot as a histogram
nbins = 10;
clc; close all
figure('color', 'white')
tiledlayout(5, 1, 'tilespacing', 'compact')


% Loop through turbines
for turbine = 1:5

    disp(turbine)
    h(turbine) = nexttile;
    ax = gca;

    % Fit a line to each trajectory
    slopes = nan(3 * num_images, 1);   
    
    for recording = 1:3
        for i = 1:num_images
            coefficients = polyfit(tracking(1).x, tracking(recording).raw(turbine).center_minimum(i,:), 1);
            slopes(i) = coefficients(1);
            clear i
        end
        clear recording
    end
    
    % Convert slopes into angles [degrees]
    wake_angles = atan(slopes);
    mean_wake_angle = mean(wake_angles, 'all', 'omitnan');

    if sign(mean_wake_angle) == 1
        color = 'green';
    elseif sign(mean_wake_angle) == -1
        color = 'red';
    else
        color = 'black';
    end

    hold on
    polarhistogram(wake_angles)
    hold off
    title(sprintf('Turbine %1.0f', turbine))
    % ylim([-30, 30])

    clear slopes wake_angles turbine clear ax
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
