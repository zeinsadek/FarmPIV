%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data path
mat_path = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks/Farm2Farm_10D_Gap';
blocks = [1,2,3];

% Figure path
figure_path = '/Users/zeinsadek/Desktop/Experiments/Farm/figures';

for b = 1:length(blocks)
    path = fullfile(mat_path, strcat("Block_", num2str(blocks(b)), "_APPENDED.mat"));
    tmp = load(path);
    tmp = tmp.combined;
    data(b) = tmp;
end

clear b tmp path


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all;
component = 'u';

ax = figure();
title(component);

% Plot combiend data
hold on
for b = 1:length(blocks)
    contourf(data(b).X, data(b).Y, data(b).(component), 100, 'linestyle', 'none')
end
hold off
colormap parula
% colormap coolwarm
% clim([0, 2.5])
% clim([1, 8])

% Plot turbines
color = 'black';
lineWidth = 1;
nacelleLength = 0.25;
for i = 1:5
    % Turbine Positions
    center = -9 + 3 * (i - 1) ;

    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', lineWidth)

    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', lineWidth)
end
hold off
axis equal
colorbar

% Seams
% Block 1
xline(1)
xline(4)
xline(7)
xline(10)

% Block 2
xline(1 + 20)
xline(4 + 20)
xline(7 + 20)
xline(10 + 20)

% Block 2
xline(1 + 40)
xline(4 + 40)
xline(7 + 40)
xline(10 + 40)

% Limits
xlim([-1, 51])
ylim([-13, 4.5])

set(gca, 'YDir','reverse')

% Ticks
yticks(-12:3:12)
xticks([0,1:3:10, 21:3:30, 41:3:50])

% Labels
fontSize = 16;
xlabel("$x / D$", "interpreter", "latex", "FontSize", fontSize)
ylabel("$y / D$", "interpreter", "latex", "FontSize", fontSize)

% Save
% file_name = strcat('SingleFarm_AllBlocks_', component, '.png');
% exportgraphics(ax, fullfile(figure_path, file_name), 'resolution', 300);
% close all