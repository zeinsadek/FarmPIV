%%% Video Generator for Time-Resolved PIV Data
% Zein Sadek
% PSU + Oldenburg

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions/Inpaint_nans/Inpaint_nans');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;

% Data paths
experiment  = 'Farm2Farm_20D_Gap';
recording   = 'Plane_3_Recording_1';

project     = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm';
data_folder = fullfile(project, 'data', experiment);
save_folder = fullfile(project, 'movies', experiment);

if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

data_path   = fullfile(data_folder, strcat(recording, '_DATA.mat'));
data        = load(data_path);
data        = data.output;

% Coordinates
X = data.X;
Y = data.Y;

% Means
U = data.U;
V = data.V;



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MOVIE (u)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num_images = 1200;
FPS        = 60;
levels     = 500;

v = VideoWriter(fullfile(save_folder, strcat(recording, '_U_MOVIE')),'MPEG-4');
v.FrameRate = FPS;
open(v)

clc; close all;
for i = 1:num_images
    progressbarText(i/num_images);
    ax = figure('Position', [300,300,200,600], 'Visible', 'off');

    % Remove Ticks
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);

    % Flip y-axis
    set(gca, 'YDir','reverse')

    % Make transparent
    % set(gcf, 'color', 'none');   
    % set(gca, 'color', 'none');

    hold on
    colormap(ax, jet)
    contourf(X, Y, inpaint_nans(double(U(:,:,i))), levels, 'linestyle', 'none')
    axis equal
    axis tight
    xlim([-120,120])
    ylim([-150,1250])
    clim([0, 8])
    c = colorbar();
    hold off

    frame = getframe(ax);
    close all
    writeVideo(v,frame);
end
close(v);













