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
experiment  = 'SingleFarm';
% recording   = 'Plane_1_Recording_1';
plane = 'Plane_1';

% Load tracking
for recording = 1:3
    tmp = load(fullfile(project_folder, 'tracking', experiment, strcat(plane, '_Recording_', num2str(recording), '_WAKE_TRACKING.mat')));
    tracking(recording) = tmp.output;
    clear tmp recording
end

% Constants
u_inf  = 7.5;
Sy     = 3;
D      = 80;

% PIV FPS
PIV_FPS = 1385;
dt = 1/PIV_FPS;
num_images = length(tracking(1).raw(1).center_gaussian);
time = linspace(0, dt * num_images, num_images);

clear data_folder data_path dt


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRA OF MIDDLE TURBINE AT DIFFERENT X LOCATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Even for this I will be averaging the spectra across the three recordings

x = tracking(1).x / D;

% Where to probe for signals
x_locations = [1,2,3,4];

% Tuning PSD parameters
window = hamming(round(num_images / 4));
noverlap = round(0.75 * length(window));
nfft = 2^nextpow2(length(window));

% Plot
figure('color', 'white')
hold on
% Loop through x-locations
for location = 1:length(x_locations)

    % Index location
    x_location = x_locations(location);
    [~, index] = min(abs(x  - x_location));
    label = sprintf('$x / D = %2.0f$', x_location);

    % Average over recordings
    PSDs = [];
    for recording = 1:3
    
        % Gather signal and remove mean
        signal = tracking(recording).raw(4).center_minimum(:, index);
        signal = signal - mean(signal, 'all', 'omitnan');
        
        % Take PSD
        [PSD, f] = pwelch(signal, window, noverlap, nfft, PIV_FPS);

        % Gather PSD
        PSDs(:, recording) = PSD;
    end

    % Average PSD
    PSD = mean(PSDs, 2, 'omitnan');
    
    % Convert frequency to strouhal
    St = (f * D * 1E-3) / u_inf;
    
    % Plot
    loglog(St, PSD, 'linewidth', 2, 'displayname', label)
    
end
hold off

xscale('log')
yscale('log')

legend('location', 'northeast', 'interpreter', 'latex')
% xlabel('Frequency [Hz]')
% ylabel('Power Spectral Density [units^2/Hz]')
xlabel('Strouhal Number')
ylabel('Power Spectral Density')
grid on


clear f index label location nfft noverlap PSD signal St window x_location x_locations


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECTRA OF ALL TURBINES AT FIXED X-LOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Where to probe for signals
x_location = 4;
[~, index] = min(abs(x  - x_location));

% Tuning PSD parameters
window = hamming(round(num_images / 4));
noverlap = round(0.5 * length(window));
nfft = 2^nextpow2(length(window));

clc; close all
figure('color', 'white')
title(sprintf('$x/ D = %2.0f$', x_location), 'Interpreter', 'latex')
hold on
% Loop through turbines
for turbine = 1:5

    % Label turbine
    label = sprintf('Turbine %1.0f', turbine);

    % Average over recordings
    PSDs = [];
    for recording = 1:3
    
        % Gather signal and remove mean
        signal = tracking(recording).raw(turbine).center_minimum(:, index);
        signal = signal - mean(signal, 'all', 'omitnan');
        
        % Take PSD
        [PSD, f] = pwelch(signal, window, noverlap, nfft, PIV_FPS);

        % Gather PSD
        PSDs(:, recording) = PSD;
    end

    % Average PSD
    PSD = mean(PSDs, 2, 'omitnan');
    
    % Convert frequency to strouhal
    St = (f * D * 1E-3) / u_inf;
    
    % Plot
    loglog(St, PSD, 'linewidth', 2, 'displayname', label)
    
end
hold off

xscale('log')
yscale('log')
xlim([1E-2, 1E1])
ylim([1E-3, 2E2])

legend('location', 'northeast', 'interpreter', 'latex')
xlabel('Strouhal Number')
ylabel('Power Spectral Density')
grid on


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPATIALLY AVERAGED SPECTRA (Across 3 Recordings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How many columns (x-locations)
columns = length(x);

% Tuning PSD parameters
window = hamming(round(num_images / 8));
noverlap = round(0.75 * length(window));
nfft = 2^nextpow2(length(window));

clc; close all
figure('color', 'white')
hold on

% Loop through turbines
for turbine = 1:5

    label = sprintf('Turbine %1.0f', turbine);
    
    PSD_sum_x = zeros(nfft/2 + 1, 1);  % Store running sum of averaged PSDs across x

    % Loop through x-locations
    for i = 1:columns
        PSD_stack_recordings = zeros(nfft/2 + 1, 3);  % Each column = PSD from 1 recording
        
        for recording = 1:3
            % Extract signal at this turbine, x-location, and recording
            signal = tracking(recording).raw(turbine).center_minimum(:, i);
            signal = signal - mean(signal, 'omitnan');

            % Compute PSD
            [PSD, f] = pwelch(signal, window, noverlap, nfft, PIV_FPS);
            PSD_stack_recordings(:, recording) = PSD;
        end
        
        % Average PSD across recordings for this x-location
        PSD_avg_recordings = mean(PSD_stack_recordings, 2, 'omitnan');

        % Accumulate for spatial average
        PSD_sum_x = PSD_sum_x + PSD_avg_recordings;
    end

    % Final spatial average across x
    PSD_avg_spatial = PSD_sum_x / columns;

    % Convert frequency to Strouhal number
    St = (f * D * 1E-3) / u_inf;

    % Plot
    loglog(St, PSD_avg_spatial, 'linewidth', 2, 'displayname', label)
end

hold off
xscale('log')
yscale('log')
xlim([1E-2, 1E1])
ylim([1E-3, 2E2])

legend('location', 'northeast', 'interpreter', 'latex')
xlabel('Strouhal Number')
ylabel('Power Spectral Density')
grid on




%% Tortuosity


dx = mean(diff(x));

figure('color', 'white')
hold on

for turbine = 1:5
    data = tracking(1).filtered(turbine).center_minimum;
    
    for i = 1:num_images
        wake_trajectory = data(i,:) / D;
    
        dy = diff(wake_trajectory);
        cord_length = sum(hypot(dx, dy));
        total_length = sum(hypot(x(end) - x(1), wake_trajectory(end) - wake_trajectory(1)));
    
        tortuosity(i) = cord_length / total_length;
    end

    plot(time, tortuosity, 'displayname', sprintf('Turbine %1.0f', turbine), 'linewidth', 2)

    fprintf('Aveage Tortuosity for Turbine %1.0f is %1.3f\n', turbine, mean(tortuosity))
end
hold off

legend()
xlim([0, max(time)])


% x = [1 2 3 4 5]; % Example x coordinates
% y = [2 4 1 3 2]; % Example y coordinates
% 
% dx = diff(x);
% dy = diff(y);
% d = hypot(dx, dy);
% path_length = sum(d);
% 
% disp(['Total path length: ', num2str(path_length)]);


%%


turbine = 1;




% Setup
num_x = length(x);          % Number of x-locations (columns)
num_rec = 3;                % Number of recordings
window = hamming(round(num_images / 4));
noverlap = round(0.75 * length(window));
nfft = 2^nextpow2(length(window));
Fs = PIV_FPS;

% Initialize PSD matrix
PSD_matrix = zeros(nfft/2 + 1, num_x);

% Loop through x-locations
for i = 1:num_x
    PSD_stack = zeros(nfft/2 + 1, num_rec);

    for r = 1:num_rec
        % Extract and detrend signal
        signal = tracking(r).raw(turbine).center_minimum(:, i);
        signal = signal - mean(signal, 'omitnan');

        % Compute PSD
        [PSD, f] = pwelch(signal, window, noverlap, nfft, Fs);
        PSD_stack(:, r) = PSD;
    end

    % Average over recordings
    PSD_matrix(:, i) = mean(PSD_stack, 2, 'omitnan');
end

% Convert frequency to Strouhal
St = (f * D * 1E-3) / u_inf;   % Same length as rows of PSD_matrix


figure('color', 'white');
imagesc(x, St, PSD_matrix); 
set(gca, 'YScale', 'log', 'XDir', 'normal');
colormap(turbo);
colorbar;

xlabel('$x / D$', 'interpreter', 'latex');
ylabel('Strouhal Number');
title(['Wake Centerline PSD – Turbine ', num2str(turbine)]);


% Optionally smooth or log-scale:
figure;
surf(x, St, log10(PSD_matrix), 'EdgeColor', 'none');
xlabel('$x/D$', 'Interpreter', 'latex');
ylabel('Strouhal Number');
zlabel('log_{10}(PSD)');
view(2); colorbar; colormap(turbo);

