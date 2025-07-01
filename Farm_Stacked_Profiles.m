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
% Plot Profiles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Linewidth and fontsize
lw = 2;
fs = 18;

%%% Block 1 = 1D - 10D
%%% Block 2 = 21D - 30D
%%% Block 3 = 41D - 50D

% Take profiles at specific locations (D)
x_locations = [1,4,7,10, 21,24,27,30, 41,44,47,50];

% Color code the profiles
colors = parula(50);

% Freestream
u_inf = 7.5;

% Which component to plot
component = 'uv';


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


    % Loop through experiments (redundant in this code)
    for e = 1:length(experiments)

        % Select experiment
        experiment = experiments{e};

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

        % Line label
        label = sprintf("$x / D = %2.0f$", x_location);

        % Select which column to plot
        [~, index] = min(abs(x - x_location));
        
        % Plot data
        plot(plot_data(:,index), y, 'color', colors(x_location, :), 'displayname', label, 'linewidth', lw)
    end
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
xlabel(x_label, 'Interpreter', 'latex', 'FontSize', fs)

% Get rid of variables (not needed, saves RAM)
clear ax block tile_counter colors e experiment h index leg location lw i x_location
clear start step stop t plot_data u_inf vis x y fs norm num_tiles x_label label component