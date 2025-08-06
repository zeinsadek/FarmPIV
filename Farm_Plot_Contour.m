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
experiments = {'SingleFarm'};
blocks = [1,2,3];

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

u_inf = 7.5;
experiment = 'SingleFarm';
component = 'uv';

% Font sizes
labelFontSize = 16;
turbineLineWidth = 2;

% Start figure
clc;
ax = figure('color', 'white', 'position', [400, 400, 800, 400]);
hold on

% Rename v to w for plots
if strcmp(component, 'v')
    component_label = 'w';
elseif strcmp(component, 'uv')
    component_label = 'uw';
elseif strcmp(component, 'vv')
    component_label = 'ww';
else
    component_label = component;
end

% What to normalize by
if ismember(component, {'uu', 'vv', 'uv'})
    norm = u_inf ^ 2;
    label = sprintf('$%s / u_\\infty^2$', component_label);
else
    norm = u_inf;
    label = sprintf('$%s / u_\\infty$', component_label);
end

% Plot combiend data
for b = 1:length(blocks)
    contourf(data.(experiment)(b).X, data.(experiment)(b).Y, imgaussfilt(data.(experiment)(b).(component) / norm, 3), 500, 'linestyle', 'none')
    clear b
end

% Which colormap
if ismember(component, {'v', 'uv'})
    colormap coolwarm
else
    colormap parula
end

% Colorbar
C = colorbar;
C.Label.String = label;
C.Label.Interpreter = 'latex';
C.Label.FontSize = labelFontSize;

% Plot turbines
color = 'black';
nacelleLength = 0.25;
for i = 1:5

    % Turbine Positions
    center = -9 + 3 * (i - 1) ;

    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', turbineLineWidth)

    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', turbineLineWidth)

    % Center line
    % yline(center, 'linestyle', '--')

    clear i
end
hold off


% Limits
axis equal
xlim([-1, 51])
ylim([-12, 4.5])

% Ticks
yticks(-12:3:12)
xticks([0, 1:3:10, 21:3:30, 41:3:50])

% Labels
xlabel("$x / D$", "interpreter", "latex", "FontSize", labelFontSize)
ylabel("$z / D$", "interpreter", "latex", "FontSize", labelFontSize)

% Orient z in the correct direction
set(gca, 'YDir','reverse')

% Clear RAM
clear C center color experiments fontSize label labelFontSize lineWidth
clear nacelleLength norm offset turbineLineWidth u_inf


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make folders
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/figures/SingleFarm/Contours';
save_folder = fullfile(figure_folder, experiment);
save_name_png = sprintf('%s_%s.png', experiment, component_label);
save_name_fig = sprintf('%s_%s.fig', experiment, component_label);

if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end

% Saving figure
clc; fprintf('Saving figure...\n')
exportgraphics(ax, fullfile(save_folder, save_name_png), 'Resolution', 300)
savefig(ax, fullfile(save_folder, save_name_fig))
close all; clc;
fprintf('Done saving!\n')


