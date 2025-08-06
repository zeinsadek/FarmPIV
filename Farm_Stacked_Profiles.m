%%% Compare profiles of a single case but have profiles be on one figure
%%% for a single quantity
% Zein Sadek

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
experiment = 'SingleFarm';
blocks = [1,2,3];

% Paths
blocks_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks');

% Store all data in structure
for b = 1:length(blocks)
    block = blocks(b);

    path = fullfile(blocks_path, experiment, sprintf("Block_%1.0f_APPENDED.mat", block));
    tmp = load(path);
    data.(experiment)(block) = tmp.combined;

end

clear e p b tmp path block blocks_path

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPROXIMATELY LOG-SPACED X-LOCATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Which component to plot
component = 'uv';

% Freestream
u_inf = 7.5;


%%% Block 1 = 1D - 10D
%%% Block 2 = 21D - 30D
%%% Block 3 = 41D - 50D

%%% ChatGPT where to take slices
% Define total span
x_min = 1;
x_max = 50;

% Create log-spaced positions
% You can adjust this for more/fewer slices
n_slices = 25; 
x_log = logspace(log10(x_min), log10(x_max), n_slices);

% Round to nearest integers
x_log_int = unique(round(x_log));

% Now keep only those values within your available data blocks
available_blocks = [1:10, 21:30, 41:50];
x_slices = intersect(x_log_int, available_blocks);

clc; disp('Suggested x locations for spanwise profiles:')
disp(x_slices)

% Take profiles at specific locations (D)
% x_locations = [1,4,7,10, 21,24,27,30, 41,44,47,50];
x_locations = x_slices;




%%% Plot where the slices are being taken
% Font sizes
labelFontSize = 16;
turbineLineWidth = 2;

% Start figure
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

    clear i
end

% Plot where slices are taken
for slice = 1:length(x_locations)
    xline(x_locations(slice))
    clear slice
end

% Plot where all log spaced profiles would be
for slice = 1:length(x_log_int)
    xline(x_log_int(slice), 'linestyle', '--', 'Alpha', 0.2)
    clear slice
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
clear C center color fontSize labelFontSize lineWidth available_blocks
clear fs x_log x_log_int x_max x_min x_slices
clear nacelleLength norm offset turbineLineWidth

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Color code the profiles
colors = parula(50);

% Linewidth and fontsize
lw = 2;
fs = 18;

% Plot
clc; close all; 
ax = figure('color', 'white');
hold on

% Flip y-axis to match coordinates
set(gca, 'YDir','reverse')

% Loop through x-locations
for location = 1:length(x_locations)

    % Where to take slice
    x_location = x_locations(location);

    % Automatically select which block
    if x_location >= 1 && x_location <= 10
        block = 1;
    elseif x_location >= 21 && x_location <= 30
        block = 2;
    elseif x_location >= 41 && x_location <= 50
        block = 3;
    else
        fprintf("Sorry buster, we dont have data at X/D = %2.0f :(\n", x_location)
        break
    end
    fprintf("Plotting data at X/D = %2.0f from block %1.0f\n", x_location, block)


    % Normalize differently based on what we are plotting
    if ismember(component, {'u', 'v'})
        norm = u_inf;
        % x_label = sprintf("$\\overline{%s} / u_{\\infty}$", component);
    else 
        norm = u_inf^2;
        % x_label = sprintf("$\\overline{%s} / u_{\\infty}^2$", component);
    end
    
    % Load data
    plot_data = data.(experiment)(block).(component) / norm;
    x = data.(experiment)(block).X(1,:);
    y = data.(experiment)(block).Y(:,1);

    % Line label
    line_label = sprintf("$x / D = %2.0f$", x_location);

    % Select which column to plot
    [~, index] = min(abs(x - x_location));
    
    % Plot data
    plot(plot_data(:,index), y, 'color', colors(x_location, :), 'displayname', line_label, 'linewidth', lw)
end


hold off

% Set y-ticks to where turbines are
yticks(-12:3:3)

% Legend
legend('location', 'northeastoutside', 'interpreter', 'latex', 'box', 'on')

% Mark where turbines are
for t = 1:5
    P = yline(-12 + 3*t, 'HandleVisibility', 'off', 'color', 'black', 'linestyle', '--', 'layer', 'bottom');
    P.Color(4) = 0.2;
end

% Axes limits
ylim([-12, 4.5])

% Axes labels
ylabel('$z / D$', 'Interpreter', 'latex', 'FontSize', fs)
xlabel(label, 'Interpreter', 'latex', 'FontSize', fs)

% Get rid of variables (not needed, saves RAM)
% clear ax block tile_counter colors e h index leg location lw i x_location
% clear start step stop t plot_data u_inf vis x y fs norm num_tiles x_label label component


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE PROFILE PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make folders
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/figures/SingleFarm/Profiles';
save_folder = fullfile(figure_folder, experiment);
save_name_png = sprintf('%s_%s_%2.0f_Slices.png', experiment, component_label, n_slices);
save_name_fig = sprintf('%s_%s_%2.0f_Slices.fig', experiment, component_label, n_slices);

if ~exist(save_folder, 'dir')
    mkdir(save_folder)
end

% Saving figure
clc; fprintf('Saving figure...\n')
exportgraphics(ax, fullfile(save_folder, save_name_png), 'Resolution', 300)
savefig(ax, fullfile(save_folder, save_name_fig))
close all; clc;
fprintf('Done saving!\n')

