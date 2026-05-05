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
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
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

% Test method by looking at single turbine
u_inf  = 8;
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
contourf(X / D, Y / D, inpaint_nans(U(:,:,1)), 20, 'linestyle', 'none')
axis equal
yticks(-12:3:3)
xticks(1:4)
clim([0.2, 1])
colorbar()
colormap(slanCM('greens'))
set(gca, 'YDir', 'reverse')
% yline(4.5)
% yline(1.5)
% yline(-1.5)
% yline(-4.5)
% yline(-7.5)
% yline(-10.5)

clear left_crop right_crop top_crop left_crop_index right_crop_index top_crop_index bottom_crop bottom_crop_index


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MOVIE BY SAVING FRAMES AS FIGURES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_folder = '/Users/zeinsadek/Downloads/EuroMech_Figures';

% Video settings
FPS = 10;
skip = 10;

% Temporary frame folder
tmp_folder = fullfile(save_folder, 'tmp_frames');
if ~exist(tmp_folder,'dir'), mkdir(tmp_folder); end

% Mac Laptop is 1440 x 900
dpi = 2 * 96;
target_width  = 300;
target_height = 1120;

% Fontsizes
labelFontSize = 10;
tickFontSize = 8;


% Generate and save frames
num_images = 2000;
c = 1;
clc;
for k = 1:skip:num_images
    progressbarText(k/num_images);

    % Update contour data only (fast)
    fig = figure('Visible', 'off', 'Color', 'white', 'units', 'inches', 'OuterPosition', [1 1 target_width/dpi target_height/dpi], 'resize', 'off');
    tiledlayout(1, 1, 'padding', 'tight');
    tile = nexttile;

    % Plot
    clc
    hold on
    contourf(X / D, Y / D, inpaint_nans(filloutliers(U(:,:,k), nan)), 20, 'linestyle', 'none')
    hold off

    axis equal
    xlim([1,4]); 
    ylim([-12, 4.5])
    set(tile, 'TickLabelInterpreter', 'latex')
    set(tile, 'FontSize', tickFontSize)

    colormap(slanCM('greens'))
    clim([0, 1])
    set(gca, 'YDir', 'reverse'); 
    yticks(-12:3:3)
    xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'FontSize', labelFontSize)
    ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'FontSize', labelFontSize)

    C = colorbar;
    C.Label.FontSize = labelFontSize;
    C.Label.String = '$u \mathbin{/} u_{\infty}$';
    C.Label.Interpreter = 'latex';
    C.TickLabelInterpreter = 'latex';
    C.FontSize = tickFontSize;
     
    % Save frame cleanly using exportgraphics
    framefile = fullfile(tmp_folder, sprintf('frame_%04d.png', c));
    exportgraphics(fig, framefile, 'Resolution', dpi);
    c = c + 1;
    close(fig)
    pause(0.1)
end



% Assemble video
vfile = fullfile(save_folder, 'SingleFarm_Plane_1_Recording_1_WakeCenter_EuroMech.mp4');
v = VideoWriter(vfile, 'MPEG-4');
v.FrameRate = FPS;

open(v)
frames = dir(fullfile(tmp_folder, '*.png'));

for k = 1:length(frames)
    img = imread(fullfile(tmp_folder, frames(k).name));
    writeVideo(v, img);
end
close(v)


% Clean up
clc;
disp('Movie generation complete.');
fprintf('Saved video: %s\n', vfile);
if exist(tmp_folder, "dir")
    rmdir(tmp_folder, "s");  % "s" = recursive delete
end
fprintf('Temporary frames folder deleted!\n\n')


