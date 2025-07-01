%% Farm2Farm Velocity spectra
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
U = double(data.U);

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
X = X(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index) / D;
Y = Y(top_crop_index:bottom_crop_index, left_crop_index:right_crop_index) / D;

x = X(1,:);
y = Y(:,1);

figure()
contourf(X, Y, U(:,:,1), 100, 'linestyle', 'none')
axis equal
yticks(-12:3:3)
colorbar()
set(gca, 'YDir', 'reverse')
yline(1.5)
yline(-1.5)
yline(-4.5)
yline(-7.5)
yline(-10.5)

clear left_crop right_crop top_crop left_crop_index right_crop_index top_crop_index bottom_crop bottom_crop_index


%% get signal at a specific point


% PIV timing
num_images = data.D;
PIV_FPS = 1385;
dt = 1/PIV_FPS;
time = linspace(0, dt * num_images, num_images);

% Where to take signal
x_position = 1;
y_position = 0;

% Indicies for position
[~, X_index] = min(abs(x - x_position));
[~, Y_index] = min(abs(y - y_position));

% Get signal and fill in nans
u_signal = squeeze(U(Y_index, X_index, :));
u_signal = fillmissing(u_signal, 'spline');
u_signal_fluctuations = u_signal - mean(u_signal, 'all', 'omitnan');

% Plot
figure()
hold on
plot(time, u_signal, 'DisplayName', 'Velocity Signal');
plot(time, u_signal_fluctuations, 'DisplayName', 'Velocity Fluctuations');
hold off
xlim([0, max(time)])
legend('location', 'NorthEast')
xlabel('Time [s]')
ylabel('u [m/s]')

clear dt X_index Y_index

%% Compute spectra

% Tuning PSD parameters
% window = hamming(256);
window = hamming(round(num_images / 12));
% noverlap = 128;
noverlap = round(0.5 * length(window));
% nfft = 512;
nfft = 2^nextpow2(length(window));


[PSD, f] = pwelch(u_signal_fluctuations, window, noverlap, nfft, PIV_FPS);

% Plot
figure;
loglog(f, PSD)
xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [units^2/Hz]')
title('Turbulence Power Spectral Density')
grid on


%% Loop through different turbines 

% Where in x to probe
x_position = 1;
[~, X_index] = min(abs(x - x_position));

% Probe right behind turbines
y_positions = -9:3:3;

% Tuning PSD parameters
window = hamming(round(num_images / 12));
noverlap = round(0.5 * length(window));
nfft = 2^nextpow2(length(window));

% Plot
clc; close all;
figure('color', 'white')
hold on
for turbine = 1:5

    % Indicies for position    
    [~, Y_index] = min(abs(y - y_positions(turbine)));
    
    % Get signal and fill in nans
    u_signal = squeeze(U(Y_index, X_index, :));
    u_signal = fillmissing(u_signal, 'spline');
    u_signal_fluctuations = u_signal - mean(u_signal, 'all', 'omitnan');

    % Compute PSD
    [PSD, f] = pwelch(u_signal_fluctuations, window, noverlap, nfft, PIV_FPS);

    % Plot
    loglog(f, PSD, 'displayname', sprintf('Turbine %1.0f', turbine))

end

loglog(f, f.^(-5/3), 'color', 'black')

hold off
xscale('log')
yscale('log')
legend()
grid on

xlabel('Frequency [Hz]')
ylabel('Power Spectral Density [units^2/Hz]')

