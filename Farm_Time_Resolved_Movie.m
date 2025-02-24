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
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/colormaps');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT PARAMETERS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
% Data paths
data_folder = '/Volumes/Zein_PIV_4/Oldenburg_Convergence/data';
recording   = 'Plane_1_Recording_3';
save_folder = '/Volumes/Zein_PIV_4/Oldenburg_Convergence/movies';

data_path   = fullfile(data_folder, strcat(recording, '_DATA.mat'));
data        = load(data_path);
data        = data.output;

% Coordinates
X = data.X.';
Y = data.Y.';

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

clc;
for i = 1:num_images
    progressbarText(i/num_images);
    ax = figure('Position', [300,300,200,600], 'Visible', 'off');
    % Remove Ticks
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);

    % Make transparent
    % set(gcf, 'color', 'none');   
    % set(gca, 'color', 'none');

    hold on
    colormap(ax, jet)
    contourf(X, Y, inpaint_nans(double(U(:,:,i))), levels, 'linestyle', 'none')
    axis equal
    axis tight
    xlim([-120,120])
    ylim([-1250,150])
    clim([0, 8])
    c = colorbar();
    % c.Label.String = 'm/s';
    hold off

    frame = getframe(ax);
    close all
    writeVideo(v,frame);
end
close(v);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MOVIE (v)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% v = VideoWriter(fullfile(save_folder, strcat(recording, '_V_MOVIE')),'MPEG-4');
% v.FrameRate = FPS;
% open(v)
% 
% for i = 1:num_images
%     progressbarText(i/num_images);
%     ax = figure('Position', [300,300,200,600], 'Visible', 'off');
%     % Remove Ticks
%     set(gca,'YTickLabel',[]);
%     set(gca,'XTickLabel',[]);
% 
%     % Make transparent
%     % set(gcf, 'color', 'none');   
%     % set(gca, 'color', 'none');
% 
%     hold on
%     colormap(ax, coolwarm)
%     contourf(X, Y, V(:,:,i), 500, 'linestyle', 'none')
%     axis equal
%     axis tight
%     xlim([-120,120])
%     ylim([-1250,150])
%     clim([-0.4, 0.4])
%     c = colorbar();
%     % c.Label.String = 'm/s';
%     hold off
% 
%     frame = getframe(ax);
%     close all
%     writeVideo(v,frame);
% end
% close(v);











