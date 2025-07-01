%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data path
experiments = {'SingleFarm', 'Farm2Farm_20D_Gap'};
blocks = [1,2,3];
% blocks = [1];

% Paths
blocks_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks');

% Store all data in structure
for b = 1:length(blocks)
    block = blocks(b);
    for e = 1:length(experiments)
        experiment = experiments{e};
    
        path = fullfile(blocks_path, experiment, sprintf("Block_%1.0f_APPENDED.mat", block));
        tmp = load(path);
        data.(experiment)(block) = tmp.combined;
    end
end

clear e p b tmp path experiment block blocks_path

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_inf = 8;

experiment = 'Farm2Farm_20D_Gap';
component = 'u';

ax = figure();
hold on
title(sprintf("%s: %s", experiment, component), 'interpreter', 'none');

if ismember(component, {'uu', 'vv', 'uv'})
    norm = u_inf ^ 2;
else
    norm = u_inf;
end

clim([0.3, 1])

% Plot combiend data
for b = 1:length(blocks)
    contourf(data.(experiment)(b).X, data.(experiment)(b).Y, imgaussfilt(data.(experiment)(b).(component) / norm, 3), 500, 'linestyle', 'none')
end

if ismember(component, {'v', 'uv'})
    colormap coolwarm
else
    colormap parula
end

% clim([0, 2.5])
% clim([2, 8])

% Plot turbines
color = 'black';
lineWidth = 3;
nacelleLength = 0.25;
for i = 1:5
    % Turbine Positions
    center = -9 + 3 * (i - 1) ;

    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', lineWidth)

    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', lineWidth)

    % Center line
    % yline(center', 'linestyle', '--')
end
hold off
axis equal
colorbar

% Seams
offset = 0;
xline(1 + offset)
xline(4 + offset)
xline(7 + offset)
xline(10 + offset)

% Limits
xlim([-1, 51])
% xlim([-1, 11])
ylim([-13, 4.5])

% Ticks
yticks(-12:3:12)
xticks([0,1 + offset:3:12 + offset])

% Labels
fontSize = 16;
xlabel("$x / D$", "interpreter", "latex", "FontSize", fontSize)
ylabel("$z / D$", "interpreter", "latex", "FontSize", fontSize)

set(gca, 'YDir','reverse')