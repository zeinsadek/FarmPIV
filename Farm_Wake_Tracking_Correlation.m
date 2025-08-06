% Compute spectra of wake meandering trajectories

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
experiments = {'SingleFarm', 'Farm2Farm_10D_Gap', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};
plane = 'Plane_1';

% Load tracking
for e = 1:length(experiments)
    experiment = experiments{e};
    for recording = 1:3
        tmp = load(fullfile(project_folder, 'tracking', experiment, strcat(plane, '_Recording_', num2str(recording), '_WAKE_TRACKING.mat')));
        tracking(recording).(experiment) = tmp.output;
        clear tmp recording e
    end
end

% Constants
u_inf  = 7.5;
Sy     = 3;
D      = 80;

% PIV FPS
PIV_FPS = 1385;
dt = 1/PIV_FPS;
num_images = length(tracking(1).SingleFarm.raw(1).center_gaussian);
time = linspace(0, dt * num_images, num_images);

clear data_folder data_path dt experiment

%% Test pearson correlation 

% X locations
x = tracking(1).SingleFarm.x / D;

% Which turbine to compare to 
reference_turbine = 1;

% Loop through turbines
clc; close all;
figure('color', 'white')
tiledlayout(1, length(experiments))

for e = 1:length(experiments)
    experiment = experiments{e};
    h(e) = nexttile;
    hold on
    for turbine = 2:5
        correlation_x = nan(1, length(x));
        for i = 1:length(x) - 1
        
            % Get turbine signals
            reference_signal = vertcat(tracking(1).(experiment).raw(reference_turbine).center_minimum(:, i), tracking(2).(experiment).raw(reference_turbine).center_minimum(:, i), tracking(3).(experiment).raw(reference_turbine).center_minimum(:, i));
            reference_signal = reference_signal - mean(reference_signal, 'all', 'omitnan');
            
            comparison_signal = vertcat(tracking(1).(experiment).raw(turbine).center_minimum(:, i), tracking(2).(experiment).raw(turbine).center_minimum(:, i), tracking(3).(experiment).raw(turbine).center_minimum(:, i));
            comparison_signal = comparison_signal - mean(comparison_signal, 'all', 'omitnan');
        
            % Compute correlation coefficient
            correlation_x(1, i) = pearson_corr(reference_signal, comparison_signal);
        
        end
    
        % Generate label
        label = sprintf('T%1.0f', turbine);
    
        plot(x, correlation_x, 'linewidth', 2, 'DisplayName', label)
    end
    hold off
    xlim([1,4])
    xlabel('$x / D$', 'interpreter', 'latex')
    legend('location', 'southwest')
    title(experiment, 'Interpreter', 'none')
end

linkaxes(h, 'xy')



%%

function r = pearson_corr(x, y)
    % PEARSON_CORR Computes the Pearson correlation coefficient between x and y
    %
    %   r = pearson_corr(x, y)
    %
    % Inputs:
    %   x - first input signal (vector)
    %   y - second input signal (vector)
    %
    % Output:
    %   r - Pearson correlation coefficient

    % Ensure x and y are column vectors
    x = x(:);
    y = y(:);
    
    % Check that x and y are the same length
    if length(x) ~= length(y)
        error('Inputs x and y must have the same length.');
    end

    % Remove NaNs (optional, comment out if not desired)
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    % Subtract means
    x_mean = mean(x);
    y_mean = mean(y);
    x_centered = x - x_mean;
    y_centered = y - y_mean;

    % Compute covariance
    covariance = mean(x_centered .* y_centered);

    % Compute standard deviations
    std_x = std(x);
    std_y = std(y);

    % Compute Pearson correlation coefficient
    r = covariance / (std_x * std_y);
end
