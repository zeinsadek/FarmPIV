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
% experiments = {'SingleFarm', 'Farm2Farm_10D_Gap', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};
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

%% Correct the Edges where planes are combined

% Each PIV Plane is 65-indicies wide

figure()
hold on
for e = 1:length(experiments)
    % for block = 1:3
    for block = 1:3
    
        % Pull in block data
        u = data.(experiments{e})(block).u;
        X = data.(experiments{e})(block).X;
        Y = data.(experiments{e})(block).Y;

        % Get x vector
        x = X(1,:);

        % Mark where seams are
        xline(x(65))
        xline(x(65 * 2))

        % Compare the ends of each plane to the start of the next
        blended_u = u;

        % Fix first edge (65 - 66)
        left = u(:, 65);
        right = u(:, 66);
        discontinuity = left - right;
        blended_u(:, 66:end) = u(:, 66:end) + discontinuity;

        clear left right discontinuity

        % Fix second edge (130 - 131)
        left = blended_u(:, 130);
        right = blended_u(:, 131);
        discontinuity = left - right;
        blended_u(:, 131:end) = blended_u(:, 131:end) + discontinuity;

        clear left right discontinuity

        % Remove crazy edge values near perimeter of PIV FOV
        blended_u(Y < -12) = nan;

        % Save fixed arrays
        cleaned.(experiments{e})(block).u = blended_u;
        cleaned.(experiments{e})(block).X = X;
        cleaned.(experiments{e})(block).Y = Y;

        % Plot
        contourf(X, Y, blended_u, 40, 'linestyle', 'none')

    end
    clear e
end
hold off

axis equal
set(gca, 'YDir', 'reverse')
colorbar()
ylim([-12, 4.5])
% clim([1, 10])


%% Compare contours

levels = 5;

figure('color', 'white')
tiledlayout(2,1)

h(1) = nexttile;
title('Original Data')
% Plot un-corrected data
hold on
for e = 1:length(experiments)
    for block = 1:3
        u = data.(experiments{e})(block).u;
        X = data.(experiments{e})(block).X;
        Y = data.(experiments{e})(block).Y;

        u(Y < -12) = nan;

        x = X(1,:);
        contourf(X, Y, u, levels, 'linestyle', 'none')
        xline(x(65))
        xline(x(65 * 2))
    end
    clear e
end
hold off

axis equal
set(gca, 'YDir', 'reverse')
colorbar()
ylim([-12, 4.5])

h(2) = nexttile;
title('Massaged Data')
% Plot corrected data
hold on
for e = 1:length(experiments)
    for block = 1:3
        u = cleaned.(experiments{e})(block).u;
        X = cleaned.(experiments{e})(block).X;
        Y = cleaned.(experiments{e})(block).Y;

        x = X(1,:);
        contourf(X, Y, u, levels, 'linestyle', 'none')
        xline(x(65))
        xline(x(65 * 2))
    end
    clear e
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
colorbar()
ylim([-12, 4.5])


% Find all axes in the current figure
% ax = findall(gcf, 'Type', 'axes');
% clims = cell2mat(get(ax, 'CLim'));
% global_clim = [min(clims(:,1)), max(clims(:,2))];
% global_clim = [0, 10];
% set(ax, 'CLim', global_clim);


%% Test cropping the edge to see how they look

levels = 50;
u_inf = 7.5;

figure('color', 'white')
tiledlayout(2,1)

h(1) = nexttile;
title('Original Data')
% Plot un-corrected data
hold on
for e = 1:length(experiments)
    for block = 1:3
        u = data.(experiments{e})(block).u;
        X = data.(experiments{e})(block).X;
        Y = data.(experiments{e})(block).Y;

        u(Y < -12) = nan;
        u(u > 0.99 * u_inf) = nan;

        x = X(1,:);
        contourf(X, Y, u, levels, 'linestyle', 'none')
        xline(x(65))
        xline(x(65 * 2))
    end
    clear e
end
hold off

axis equal
set(gca, 'YDir', 'reverse')
colorbar()
ylim([-12, 4.5])

h(2) = nexttile;
title('Massaged Data')
% Plot corrected data
hold on
for e = 1:length(experiments)
    for block = 1:3
        u = cleaned.(experiments{e})(block).u;
        X = cleaned.(experiments{e})(block).X;
        Y = cleaned.(experiments{e})(block).Y;

        u(u > 0.99 * u_inf) = nan;

        x = X(1,:);
        contourf(X, Y, u, levels, 'linestyle', 'none')
        xline(x(65))
        xline(x(65 * 2))
    end
    clear e
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
colorbar()
ylim([-12, 4.5])


%% Loop and correct all components

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Save paths
save_path = '/Users/zeinsadek/Downloads';

% Each PIV Plane is 65-indicies wide
for e = 1:length(experiments)
    for block = 1:3
        for c = 1:length(components)

            component = components{c};
            % Pull in block data
            u = data.(experiments{e})(block).(component);
            X = data.(experiments{e})(block).X;
            Y = data.(experiments{e})(block).Y;
    
            % Get x vector
            x = X(1,:);

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
            cleaned.(experiments{e})(block).(component) = blended;
            cleaned.(experiments{e})(block).X = X;
            cleaned.(experiments{e})(block).Y = Y;
        end
    end
    
    % Save to matfile
    save(fullfile(save_path, sprintf('%s_Massaged_APPENDED.mat', experiments{e})), 'cleaned')

    clear e cleaned u X Y
end



