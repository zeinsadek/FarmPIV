%%% Compare profiles between different cases
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
experiments = {'SingleFarm', 'Farm2Farm_10D_Gap', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};
% experiments = {'SingleFarm', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};
blocks = [1,2,3];
% blocks = 1;

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

% Load incident inflows
inflows = load('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/matfiles/SingleFarm_InformedInflows.mat');
inflows = inflows.output;

clear e p b tmp path experiment block blocks_path

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot Profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Colors
colors.SingleFarm = 'black';
colors.Farm2Farm_10D_Gap = 'red';
colors.Farm2Farm_20D_Gap = 'green';
colors.Farm2Farm_40D_Gap = 'blue';

% Linewidth and fontsize
lw = 2;
fs = 16;

%%% Block 1 = 1D - 10D
%%% Block 2 = 21D - 30D
%%% Block 3 = 41D - 50D

% Freestream
% u_inf = 8;

% Take profiles at specific locations (D)
start = 1;
stop = 10;
step = 1;

% Which component to plot
component = 'u';

% Automatically select which block
if start == 1
    block = 1;
elseif start == 21
    block = 2;
elseif start == 41
    block = 3;
end
clc; fprintf("Plotting data from block %1.0f\n", block)


% Plot
close all; 
ax = figure('color', 'white');
num_tiles = (stop - start + 1) / step;
tiledlayout(1, num_tiles)

% Counter for tiles
tile_counter = 1;

% Loop through x-locations
for location = start:step:stop

    % Start plot
    h(tile_counter) = nexttile;
    hold on

    % Flip y-axis to match coordinates
    set(gca, 'YDir','reverse')

    % Loop through experiments
    for e = 1:length(experiments)

        % Select experiment
        experiment = experiments{e};

        % Load incident inflow
        u_inf = inflows.(experiment).u_incident;

        % Normalize differently based on what we are plotting
        if ismember(component, {'u', 'v'})
            norm = u_inf;
            x_label = sprintf("$\\overline{%s} / u_{\\infty}$", component);
        else 
            norm = u_inf^2;
            x_label = sprintf("$\\overline{%s} / u_{\\infty}^2$", component);
        end
        
        % Load data
        plot_data = data.(experiment)(block).(component) / norm;
        x = data.(experiment)(block).X(1,:);
        y = data.(experiment)(block).Y(:,1);

        % Select which column to plot
        [~, index] = min(abs(x - location));
        
        % Only make legend work for first plot
        if tile_counter == 1
            vis = 'on';
        else 
            vis = 'off';
        end

        % Plot data
        plot(plot_data(:,index), y, 'color', colors.(experiment), 'displayname', experiment, 'linewidth', lw, 'HandleVisibility', vis)
    end
    hold off

    % Set y-ticks to where turbines are
    yticks(-12:3:3)

    % Title which x-location each plot is taken at
    title(sprintf('$x / D = %1.0f$', location), 'interpreter', 'latex')
    
    % Make a legend for only the first plot
    if tile_counter == 1
        leg = legend('Orientation', 'Horizontal', 'Interpreter', 'none');
        leg.Layout.Tile = 'north';
        leg.Box = 'off';
        ylabel('$z / D$', 'Interpreter', 'latex', 'fontsize', fs)
    end

    % Mark where turbines are
    % for t = 1:5
    %     yline(-12 + 3*t, 'HandleVisibility', 'off')
    % end

    % Add x-label for the middle plot
    if tile_counter == round(num_tiles/2)
        xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fs)
    end

    % Hide y-axis for plots except the first one
    if tile_counter ~= 1
        h(tile_counter).YAxis.Visible = 'off';   % remove y-axis
    end

    % Increment tile counter
    tile_counter = tile_counter + 1;
end

% Make all tiles have the same x and y limits
linkaxes(h, 'xy')

% Set unified axis limits
if ismember(component, {'u', 'uu', 'vv'})
    % Set positive-definite x-limits for all linked axes
    for i = 1:length(h)
        xlim(h(i), [0, max(xlim(h(i)))])  % clip lower bound only
        ylim(h(i), [-12, 4.5])
    end
else
    % Still apply consistent y-limits
    for i = 1:length(h)
        ylim(h(i), [-12, 4.5])
    end
end


% Get rid of variables (not needed, saves RAM)
clear ax block tile_counter colors e experiment h index leg location lw i 
clear start step stop t plot_data u_inf vis x y fs norm num_tiles x_label