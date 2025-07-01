%% Find adjusted inflow for Farm2Farm cases using single farm wake as reference

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
% GET PROFILE AT SPECIFIED LOCATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_locations = [10, 20, 40];
save_names = {'Farm2Farm_10D_Gap', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};

clc; close all
% Where to take profile
for s = 1:length(x_locations)
    x_location = x_locations(s);
    
    if x_location >= 0 && x_location <= 11
        block = 1;
    elseif x_location >= 20 && x_location <= 31
        block = 2;
    elseif x_location >= 40 && x_location <= 51
        block = 3;
    end
    
    % Get correct data
    u = data.SingleFarm(block).u;
    X = data.SingleFarm(block).X;
    Y = data.SingleFarm(block).Y;
    
    x = X(1,:);
    y = Y(:,1);
    
    % Where to take slice
    [~, X_index] = min(abs(x - x_location));
    fprintf("%s: Want a slice at %2.0fD, taking slice at %2.1fD\n", save_names{s}, x_location, x(X_index))
    
    % Get velocity profile
    u_profile = u(:, X_index);
    
    % Aveage within the bounds of the farm
    top_bound = -9;
    bottom_bound = 4.5;
    [~, top_bound_index] = min(abs(y - top_bound));
    [~, bottom_bound_index] = min(abs(y - bottom_bound));
    
    % Get average velocity within farm
    u_incident = mean(u_profile(top_bound_index:bottom_bound_index), 'all', 'omitnan');
    fprintf("Incident velocity is %2.2f\n\n", u_incident)
    
    %%% Plot contour to check 
    figure('color', 'white')
    tiledlayout(1,2)
    sgtitle(save_names{s}, 'interpreter', 'none')

    nexttile
    hold on
    contourf(X, Y, u, 100, 'LineStyle', 'none')
    xline(x_location, 'color', 'black', 'linestyle', '--')
    xline(x(X_index), 'color', 'red', 'LineWidth', 2)
    set(gca, 'YDir', 'reverse')
    axis equal
    xlabel('X / D')
    ylabel('Z / D')
    xlim([min(x) - 2, max(x) + 1])
    ylim([-12, 4.5])
    yticks(-12:3:3)
    
    % Mark where turbines are
    for t = -9:3:3
        yline(t, 'linestyle', '--')
    end
    
    %%% Plot profile
    nexttile
    hold on
    plot(u_profile, y)
    plot(u_profile(top_bound_index:bottom_bound_index), y(top_bound_index:bottom_bound_index), 'color', 'red')
    
    % Plot incident velocity
    xline(u_incident, 'color', 'red', 'LineWidth', 2)
    
    hold off
    set(gca, 'YDir', 'reverse')
    ylim([-12, 4.5])
    xlabel('u [m/s]')
    ylabel('Z / D')
    xlim([6, 8.5])
    yticks(-12:3:3)
    
    % Mark where turbines are
    for t = -9:3:3
        yline(t, 'linestyle', '--')
    end

    % Save to output array
    output.(save_names{s}).u_incident = u_incident;
    output.(save_names{s}).profile = u_profile;
end

% Add single farm case
output.SingleFarm.u_incident = 8;

clear block bottom_bound bottom_bound_index closest_x s t top_bound top_bound_index u u_incident u_profile
clear x X y Y X_index x_location x_locations 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_dir = "/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/matfiles";
save_name = 'SingleFarm_InformedInflows.mat';

save(fullfile(save_dir, save_name), 'output');
clc; fprintf("Matfile saved!\n")
