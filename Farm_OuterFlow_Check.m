%%% Check outer region of flow and that it is consistant across all the
%%% different cases we have
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
blocks = [1, 2, 3];

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
% LOAD INFLOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

inflow = load('/Users/zeinsadek/Downloads/Plane_9_COMBINED.mat');
inflow = inflow.output;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE A MASSAGED VERSION OF THE DATA AT THE SEAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Loop through cases
for e = 1:length(experiments)
    experiment = experiments{e};

    % Loop through blocks
    for block = 1:3
        for c = 1:length(components)
            component = components{c};
    
            % Pull in block data
            u = data.(experiment)(block).(component);
            X = data.(experiment)(block).X;
            Y = data.(experiment)(block).Y;
    
            % Compare the ends of each plane to the start of the next
            blended = u;
    
            % Fix first edge (65 - 66)
            left = u(:, 65);
            right = u(:, 66);
            discontinuity = left - right;
            blended(:, 66:end) = u(:, 66:end) + discontinuity;
    
            clear left right discontinuity
    
            % Fix second edge (130 - 131)
            left = blended(:, 130);
            right = blended(:, 131);
            discontinuity = left - right;
            blended(:, 131:end) = blended(:, 131:end) + discontinuity;
    
            clear left right discontinuity
    
            % Remove crazy edge values near perimeter of PIV FOV
            blended(Y < -12) = nan;
    
            % Save fixed arrays
            cleaned.(experiment)(block).(component) = blended;
            cleaned.(experiment)(block).X = X;
            cleaned.(experiment)(block).Y = Y;
        end
    end
end

clear components block c component u X Y blended left right discontinuity 
clear massage_outer_top massage_outer_bottom e experiment


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT OUTER EDGE OF FLOW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Static block for outer flow
outer_top = -3;
outer_bottom = -2.75;

% Static block for inner core
inner_top = 3;
inner_bottom = 14;

% Colors for different cases
colors.SingleFarm = 'black';
colors.Farm2Farm_10D_Gap = 'red';
colors.Farm2Farm_20D_Gap = 'green';
colors.Farm2Farm_40D_Gap = 'blue';


figure('color', 'white')
title('Massaged')
hold on
for e = 1:length(experiments)
    experiment = experiments{e};
    for block = 1:3
    
        % Load data
        u = cleaned.(experiment)(block).u;
        Y = data.(experiment)(block).Y - (-9);
        X = data.(experiment)(block).X;
        y = Y(:,1);
        x = X(1,:);
    
        % Mask inner region
        [~, inner_top_idx] = min(abs(y - inner_top));
        [~, inner_bottom_idx] = min(abs(y - inner_bottom));  
        u_inner = u(inner_top_idx:inner_bottom_idx, :);
    
        % Mask outer region
        [~, outer_top_idx] = min(abs(y - outer_top));
        [~, outer_bottom_idx] = min(abs(y - outer_bottom));  
        u_outer = u(outer_top_idx:outer_bottom_idx, :);
   
        
        % Average inner/outer regions 
        inner_mean(block).u = mean(u_inner, 1, 'omitnan');
        outer_mean(block).u = mean(u_outer, 1, 'omitnan');
    

        % Get center line velocity of middle turbine
        % [~, center_idx] = min(abs(y - Sy * 3));
        % middle_centerline(block).u = u(center_idx, :);

        % plot(x, inner_mean(block).u)
        plot(x, outer_mean(block).u, 'color', colors.(experiment), 'linewidth', 2)
    
    
    end
end

% Inflow data
yline(mean(inflow.u, 'all', 'omitnan'), 'linewidth', 2)

hold off
ylim([0, 9])
xlim([0, 50])
xlabel('x / D')
ylabel('Outer Velocity [m/s]')



%% Test massaging for 20D and 40D gaps

experiment = 'SingleFarm';

figure('color', 'white')
tiledlayout(2,1)
sgtitle(experiment, 'interpreter', 'none')

h(1) = nexttile();
title('Raw')
hold on
for block = 1:3
    contourf(data.(experiment)(block).X, data.(experiment)(block).Y, data.(experiment)(block).u, 10, 'linestyle', 'none')
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
colorbar()
clim([2, 8])


h(2) = nexttile();
title('massaged')
hold on
for block = 1:3
    contourf(data.(experiment)(block).X, data.(experiment)(block).Y, cleaned.(experiment)(block).u, 10, 'linestyle', 'none')
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
colorbar()
clim([2, 8])

linkaxes(h, 'xy')