% Attempt to show the far wake can be scaled: using mixing layer approach
% Zein Sadek

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Documents/MATLAB/MatlabFunctions/Inpaint_nans')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data path
experiment = {'SingleFarm'};
blocks = [1,2,3];

% Turbine/Farm dimensions
D = 80;
Sy = 3;

% Paths
blocks_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks');

% Store all data in structure
for b = 1:length(blocks)
    block = blocks(b);

    path = fullfile(blocks_path, experiment, sprintf("Block_%1.0f_APPENDED.mat", block));
    tmp = load(path);
    data(block) = tmp.combined;

end

clear e p b tmp path experiment block blocks_path


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIX OUTER FLOW REGION IN RAW WAKE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Spacing between blocks
block_spacing = 20;

% Try to detect bad points in outer region
massage_outer_top = -4;
massage_outer_bottom = -3.5;

% How much to delete (in indicies)
mask_width = 14;

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Where turbines are   
turbine_positions = -9:3:3;


for block = 1:3
    for c = 1:length(components)
        component = components{c};

        % Load coordinates
        X = data(block).X;
        x = X(1,:);
        Y = data(block).Y - turbine_positions(1);
        y = Y(:,1);
    
        % Load data
        raw_wake = data(block).u;
    
        [~, outer_top_idx] = min(abs(y - massage_outer_top));
        [~, outer_bottom_idx] = min(abs(y - massage_outer_bottom));
     
        % Mask upper left corner of each plane
        corner_xs = [1,4,7] + (block - 1) * block_spacing;
    
        % Mask bad corners
        for j = 1:3
            [~, tmp_x_idx] = min(abs(x - corner_xs(j)));
            raw_wake(outer_top_idx:outer_bottom_idx, tmp_x_idx:(tmp_x_idx + mask_width - 1)) = nan;
        end
    
        % Fill in bad spots
        % raw_wake(outer_top_idx:outer_bottom_idx, :) = fillmissing(raw_wake(outer_top_idx:outer_bottom_idx, :), 'next', 2);

        % Try using inpaint_nans
        raw_wake = inpaint_nans(double(raw_wake));
    
        % Save fixed arrays
        outer_fixed(block).(component) = raw_wake;
        outer_fixed(block).X = X;
        outer_fixed(block).Y = Y;
    end
end


% Plot to check
figure('color', 'white')
hold on
for block = 1:3
    contourf(outer_fixed(block).X, outer_fixed(block).Y, outer_fixed(block).u, 100, 'linestyle', 'none')
end
hold off

axis equal
set(gca, 'YDir', 'reverse')
xlabel('u [m/s]')
ylabel('z / D')
ylim([-4.5, 15])

clear colors massage_outer_bottom massage_outer_top block block_spacing c component corner_xs j mask_width
clear outer_bottom_idx outer_top_idx raw_wake tmp_x_idx x X y Y components


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE A MASSAGED VERSION OF THE DATA AT THE SEAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Loop through blocks
for block = 1:3
    for c = 1:length(components)
        component = components{c};

        % Pull in block data
        u = data(block).(component);
        % u = outer_fixed(block).(component);
        X = data(block).X;
        Y = data(block).Y;

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
        cleaned(block).(component) = blended;
        cleaned(block).X = X;
        cleaned(block).Y = Y;
    end
end


% Plot to check
figure('color', 'white')
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y, cleaned(block).u, 100, 'linestyle', 'none')
end
hold off

axis equal
set(gca, 'YDir', 'reverse')
xlabel('u [m/s]')
ylabel('z / D')

clear components block c component u X Y blended left right discontinuity massage_outer_top massage_outer_bottom


%% plot profiles to try and fix fucked up edge

skip = 10;

figure('color', 'white')
hold on
for block = 1:3
    for i = 1:skip:190
        plot(data(block).u(:,i), data(block).Y(:,1))
    end
end
hold on
set(gca, 'YDir', 'reverse')
yline(-12)

%% Try fixing after massaging

% for block = 1:3
% 
%         % Load coordinates
%         X = cleaned(block).X;
%         Y = cleaned(block).Y;
%         raw_wake = cleaned(block).u;
%         x = X(1,:);
%         y = Y(:,1);
% 
% 
%         % Try to detect bad points in outer region
%         massage_outer_top = -13;
%         massage_outer_bottom = -12;
% 
%         [~, outer_top_idx] = min(abs(y - massage_outer_top));
%         [~, outer_bottom_idx] = min(abs(y - massage_outer_bottom));
% 
%         raw_wake(outer_top_idx:outer_bottom_idx, :) = nan;
% 
%         raw_wake = fillmissing(raw_wake, 'next', 1);
% 
%         % Fill in bad spots
%         % raw_wake(outer_top_idx:outer_bottom_idx, :) = fillmissing(raw_wake(outer_top_idx:outer_bottom_idx, :), 'next', 1);
% 
% 
%         % Save fixed arrays
%         test(block).u = raw_wake;
%         test(block).X = X;
%         test(block).Y = Y;
% 
% end
% 
% 
% % Plot to check
% figure('color', 'white')
% hold on
% for block = 1:3
%     contourf(test(block).X, test(block).Y, test(block).u, 100, 'linestyle', 'none')
% end
% hold off
% yline(massage_outer_top)
% yline(massage_outer_bottom)
% axis equal
% set(gca, 'YDir', 'reverse')
% xlabel('u [m/s]')
% ylabel('z / D')
% % ylim([-4.5, 15])

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOOK AT PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Near, middle, and far regions
colors(1).c = 'red';
colors(2).c = 'green';
colors(3).c = 'blue';

block_spacing = 20;

% Plot profiles
figure('color', 'white')
hold on
for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(1);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_outer = outer_fixed(block).u;
    massaged_wake = cleaned(block).u;

    % Plot
    contourf(X, Y, massaged_wake, 100, 'linestyle', 'none')

end
hold off

axis equal
set(gca, 'YDir', 'reverse')
xlabel('u [m/s]')
ylabel('z / D')
ylim([-4.5, 15])



% Plot profiles
figure('color', 'white')
hold on
for block = 1:3

    % Load coordinates
    Y = data(block).Y - turbine_positions(1);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_outer = outer_fixed(block).u;
    massaged_wake = cleaned(block).u;

    % Plot
    step = 10;
    for i = 1:step:size(raw_wake,2)
        plot(massaged_wake(:,i), y, 'color', colors(block).c)
    end
    
end
hold off
set(gca, 'YDir', 'reverse')

clear block block_spacing i massaged_outer massaged_wake raw_wake step x X y Y


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE AVG INNER AND OUTER VELOCITIES
% BLOCK AVERAGE INNER AND OUTER REGIONS (U1 & U2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Can use bigger bounds for when we can fix the outer flow
% outer_top = -4;
% outer_bottom = -3;

% Static block for outer flow
outer_top = -3;
outer_bottom = -2.75;

% Static block for inner core
inner_top = 3;
inner_bottom = 14;


figure('color', 'white')
hold on
for block = 1:3

    % Load data
    u = cleaned(block).u;
    Y = data(block).Y - turbine_positions(1);
    X = data(block).X;
    y = Y(:,1);
    x = X(1,:);

    % Mask inner region
    [~, inner_top_idx] = min(abs(y - inner_top));
    [~, inner_bottom_idx] = min(abs(y - inner_bottom));  
    u_inner = u(inner_top_idx:inner_bottom_idx, :);
    contourf(X(inner_top_idx:inner_bottom_idx, :), Y(inner_top_idx:inner_bottom_idx, :), u_inner, 20, 'LineStyle', 'none')

    % Mask outer region
    [~, outer_top_idx] = min(abs(y - outer_top));
    [~, outer_bottom_idx] = min(abs(y - outer_bottom));  
    u_outer = u(outer_top_idx:outer_bottom_idx, :);
    contourf(X(outer_top_idx:outer_bottom_idx, :), Y(outer_top_idx:outer_bottom_idx, :), u_outer, 20, 'LineStyle', 'none')

    
    % Average inner/outer regions 
    inner_mean(block).u = mean(u_inner, 1, 'omitnan');
    outer_mean(block).u = mean(u_outer, 1, 'omitnan');

    % Min/max of inner/outer regions
    inner_min(block).u = min(u_inner, [], 1, 'omitnan');
    outer_max(block).u = max(u_outer, [], 1, 'omitnan');
    xs(block).x = x;

    % Get center line velocity of middle turbine
    [~, center_idx] = min(abs(y - Sy * 3));
    middle_centerline(block).u = u(center_idx, :);


end
hold off
axis equal
set(gca, 'YDir', 'reverse')


lw = 2;
figure('color', 'white')
tile = tiledlayout(1,2);

h(1) = nexttile;
hold on
for block = 1:3

    % Means
    plot(xs(block).x, inner_mean(block).u, 'color', 'red', 'linewidth', lw, 'HandleVisibility', 'off')
    % plot(xs(block).x, outer_mean(block).u, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off')

    % Min/max
    plot(xs(block).x, inner_min(block).u, 'color', 'green', 'linestyle', '-', 'linewidth', lw, 'HandleVisibility', 'off')
    % plot(xs(block).x, outer_max(block).u, 'color', 'black', 'linestyle', '--', 'linewidth', lw, 'HandleVisibility', 'off')

    % Centerline of middle farm
    plot(xs(block).x, middle_centerline(block).u, 'color', 'blue', 'HandleVisibility', 'off', 'linewidth', lw)

end

% Legend (colors)
plot(nan, nan, 'linewidth', lw, 'color', 'red', 'DisplayName', 'Mean')
plot(nan, nan, 'linewidth', lw, 'color', 'green', 'DisplayName', 'Min')
plot(nan, nan, 'linewidth', lw, 'color', 'blue', 'DisplayName', 'Centerline')



hold off
legend('box', 'off', 'location', 'southeast')
ylabel('$u_1, \,\, u_2$ [m/s]', 'interpreter', 'latex')
xlabel('$x / D$', 'interpreter', 'latex')

h(2) = nexttile;
hold on
for block = 1:3

    % Means
    plot(xs(block).x, outer_mean(block).u - inner_mean(block).u, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off')

    % Min/max
    plot(xs(block).x, outer_max(block).u - inner_min(block).u, 'color', 'black', 'linestyle', '--', 'linewidth', lw, 'HandleVisibility', 'off')

end

% Legend
plot(nan, nan, 'linewidth', lw, 'color', 'black', 'linestyle', '-', 'DisplayName', 'Average')
plot(nan, nan, 'linewidth', lw, 'color', 'black', 'linestyle', '--', 'DisplayName', 'Min/Max')

hold off
legend('box', 'off', 'location', 'northeast')
ylabel('$\Delta u = u_1 - u_2$', 'interpreter', 'latex')
xlabel('$x / D$', 'interpreter', 'latex')

clear block h inner_bottom inner_bottom_idx inner_top inner_top_idx lw outer_bottom outer_bottom_idx outer_top outer_top_idx tile 
clear u x xs X y Y u_inner u_outer center_idx


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FARM HALF WIDTH (Mean velocity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turbine = 1;
half_width_factor = 0.5;

clc; close all
figure('color', 'white')
title('Block Mean Velocity')
hold on
for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;
    massaged_wake_uncropped = cleaned(block).u;

    % Crop inner side
    raw_wake(Y > 3) = nan;
    massaged_wake(Y > 3) = nan;
    massaged_wake_uncropped(Y > 3) = nan;

    % Get centerline velocities (WHICH U1 AND U2 WE USE)
    raw_center_velocity = inner_mean(block).u;
    massaged_center_velocity = raw_center_velocity;




    % Find half-width velocity
    u_inf = outer_mean(block).u;
    raw_center_half_velocity = half_width_factor * (raw_center_velocity + u_inf);
    massaged_center_half_velocity = half_width_factor * (massaged_center_velocity + u_inf);

    % Inside half-width
    raw_wake(raw_wake > raw_center_half_velocity) = nan;
    massaged_wake(massaged_wake > massaged_center_half_velocity) = nan;

    % Detect top and bottom half-widths
    raw_isnan = isnan(raw_wake);
    massaged_isnan = isnan(massaged_wake);



    %%% Top (bottom unflipped)
    raw_isnan_bottom = raw_isnan;
    massaged_isnan_bottom = massaged_isnan;

    % Vertical difference to find edge
    raw_isnan_bottom_diff = abs(diff(raw_isnan_bottom, 1, 1));
    massaged_isnan_bottom_diff = abs(diff(massaged_isnan_bottom, 1, 1));

    % Find half-width indicies
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);

    % Get half-width
    raw_bottom_edge = Y(raw_bottom_edge_index);
    massaged_bottom_edge = Y(massaged_bottom_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y(1);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y(1);

    % Plot
    contourf(X, Y, massaged_wake_uncropped, 'linestyle', 'none')
    plot(x, sgolayfilt(massaged_bottom_edge, 3, 5), 'color', 'black', 'linewidth', 2)
    colorbar()


    % Save values
    halfwidth_mean.x(block, :) = x;
    halfwidth_mean.massaged(block, :) = massaged_bottom_edge;
    halfwidth_mean.raw(block, :) = raw_bottom_edge;
    
end


hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-3, 12])
yline(0, 'linestyle', '--')
xlabel('$x / D$', 'interpreter', 'latex')
ylabel('$z / D$', 'interpreter', 'latex')

clear block half_width_factor massaged_bottom_edge massaged_bottom_edge_index massaged_center_half_velocity
clear massaged_center_velocity massaged_isnan massaged_isnan_bottom massaged_isnan_bottom_diff
clear massaged_wake massaged_wake_uncropped raw_bottom_edge raw_bottom_edge_index raw_center_half_velocity
clear raw_center_velocity raw_isnan raw_isnan_bottom raw_isnan_bottom_diff raw_wake turbine turbine_colors
clear x_fit y_fit x halfwidth_x X y halfwidth_massaged Y


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FARM HALF WIDTH (Min/Max velocity)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turbine = 1;
half_width_factor = 0.5;

clc; close all
figure('color', 'white')
title('Min/Max Velocities')
hold on
for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;
    massaged_wake_uncropped = cleaned(block).u;

    % Crop inner side
    raw_wake(Y > 3) = nan;
    massaged_wake(Y > 3) = nan;
    massaged_wake_uncropped(Y > 3) = nan;

    % Get centerline velocities (WHICH U1 AND U2 WE USE)
    raw_center_velocity = inner_min(block).u;
    massaged_center_velocity = raw_center_velocity;

    % Find half-width velocity
    u_inf = outer_max(block).u;



    raw_center_half_velocity = half_width_factor * (raw_center_velocity + u_inf);
    massaged_center_half_velocity = half_width_factor * (massaged_center_velocity + u_inf);

    % Inside half-width
    raw_wake(raw_wake > raw_center_half_velocity) = nan;
    massaged_wake(massaged_wake > massaged_center_half_velocity) = nan;

    % Detect top and bottom half-widths
    raw_isnan = isnan(raw_wake);
    massaged_isnan = isnan(massaged_wake);



    %%% Top (bottom unflipped)
    raw_isnan_bottom = raw_isnan;
    massaged_isnan_bottom = massaged_isnan;

    % Vertical difference to find edge
    raw_isnan_bottom_diff = abs(diff(raw_isnan_bottom, 1, 1));
    massaged_isnan_bottom_diff = abs(diff(massaged_isnan_bottom, 1, 1));

    % Find half-width indicies
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);

    % Get half-width
    raw_bottom_edge = Y(raw_bottom_edge_index);
    massaged_bottom_edge = Y(massaged_bottom_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y(1);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y(1);

    % Plot
    contourf(X, Y, massaged_wake_uncropped, 'linestyle', 'none')
    plot(x, sgolayfilt(massaged_bottom_edge, 3, 5), 'color', 'black', 'linewidth', 2)
    colorbar()


    % Save values
    halfwidth_minMax.x(block, :) = x;
    halfwidth_minMax.massaged(block, :) = massaged_bottom_edge;
    halfwidth_minMax.raw(block, :) = raw_bottom_edge;
    
end


hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-3, 12])
yline(0, 'linestyle', '--')
xlabel('$x / D$', 'interpreter', 'latex')
ylabel('$z / D$', 'interpreter', 'latex')

clear block half_width_factor massaged_bottom_edge massaged_bottom_edge_index massaged_center_half_velocity
clear massaged_center_velocity massaged_isnan massaged_isnan_bottom massaged_isnan_bottom_diff
clear massaged_wake massaged_wake_uncropped raw_bottom_edge raw_bottom_edge_index raw_center_half_velocity
clear raw_center_velocity raw_isnan raw_isnan_bottom raw_isnan_bottom_diff raw_wake turbine turbine_colors
clear x_fit y_fit x halfwidth_x X y halfwidth_massaged Y


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FARM HALF WIDTH (middle centerline)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turbine = 1;
half_width_factor = 0.5;

clc; close all
figure('color', 'white')
title('Middle Centerline Velocities')
hold on
for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;
    massaged_wake_uncropped = cleaned(block).u;

    % Crop inner side
    raw_wake(Y > 3) = nan;
    massaged_wake(Y > 3) = nan;
    massaged_wake_uncropped(Y > 3) = nan;

    % Get centerline velocities (WHICH U1 AND U2 WE USE)
    raw_center_velocity = middle_centerline(block).u;
    massaged_center_velocity = raw_center_velocity;

    % Find half-width velocity
    u_inf = outer_max(block).u;



    raw_center_half_velocity = half_width_factor * (raw_center_velocity + u_inf);
    massaged_center_half_velocity = half_width_factor * (massaged_center_velocity + u_inf);

    % Inside half-width
    raw_wake(raw_wake > raw_center_half_velocity) = nan;
    massaged_wake(massaged_wake > massaged_center_half_velocity) = nan;

    % Detect top and bottom half-widths
    raw_isnan = isnan(raw_wake);
    massaged_isnan = isnan(massaged_wake);



    %%% Top (bottom unflipped)
    raw_isnan_bottom = raw_isnan;
    massaged_isnan_bottom = massaged_isnan;

    % Vertical difference to find edge
    raw_isnan_bottom_diff = abs(diff(raw_isnan_bottom, 1, 1));
    massaged_isnan_bottom_diff = abs(diff(massaged_isnan_bottom, 1, 1));

    % Find half-width indicies
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);

    % Get half-width
    raw_bottom_edge = Y(raw_bottom_edge_index);
    massaged_bottom_edge = Y(massaged_bottom_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y(1);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y(1);

    % Plot
    contourf(X, Y, massaged_wake_uncropped, 'linestyle', 'none')
    plot(x, sgolayfilt(massaged_bottom_edge, 3, 5), 'color', 'black', 'linewidth', 2)
    colorbar()


    % Save values
    halfwidth_middle.x(block, :) = x;
    halfwidth_middle.massaged(block, :) = massaged_bottom_edge;
    halfwidth_middle.raw(block, :) = raw_bottom_edge;
    
end


hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-3, 12])
yline(0, 'linestyle', '--')
xlabel('$x / D$', 'interpreter', 'latex')
ylabel('$z / D$', 'interpreter', 'latex')

clear block half_width_factor massaged_bottom_edge massaged_bottom_edge_index massaged_center_half_velocity
clear massaged_center_velocity massaged_isnan massaged_isnan_bottom massaged_isnan_bottom_diff
clear massaged_wake massaged_wake_uncropped raw_bottom_edge raw_bottom_edge_index raw_center_half_velocity
clear raw_center_velocity raw_isnan raw_isnan_bottom raw_isnan_bottom_diff raw_wake turbine turbine_colors
clear x_fit y_fit x halfwidth_x X y halfwidth_massaged Y




%% Plot velocity deficit times the halfwidth to find where the profiles are
% roughly self-similar

lw = 2;

figure('color', 'white')
title('$\Delta U \cdot \delta_{1/2}$', 'interpreter', 'latex', 'fontsize', 24)
hold on
for block = 1:3

    % Mean velocities
    U1_mean = outer_mean(block).u;
    U2_mean = inner_mean(block).u;

    % Min/max velocities
    U1_max = outer_max(block).u;
    U2_min = inner_min(block).u;

    % Centerline velocity
    U2_middle = middle_centerline(block).u;

    % Deficit
    Delta_U_mean = U1_mean - U2_mean;
    Delta_U_minMax = U1_max - U2_mean;
    Delta_U_middle = U1_mean - U2_middle;

    % Halfwidths
    x_tmp = halfwidth_mean.x(block,:);

    delta_half_mean = halfwidth_mean.massaged(block, :);
    delta_half_minMax = halfwidth_minMax.massaged(block, :);
    delta_half_middle = halfwidth_middle.massaged(block, :);


    % Plot 
    plot(x_tmp, Delta_U_mean .* delta_half_mean, 'linewidth', lw, 'color', 'red', 'HandleVisibility', 'off')
    plot(x_tmp, Delta_U_minMax .* delta_half_minMax, 'linewidth', lw, 'color', 'green', 'HandleVisibility', 'off')
    plot(x_tmp, Delta_U_middle .* delta_half_middle, 'linewidth', lw, 'color', 'blue', 'HandleVisibility', 'off')

    % Plot horizontal line that is mean of far reigion
    if block == 3
        yline(mean(Delta_U_mean .* delta_half_mean), 'color', 'red', 'HandleVisibility', 'off')
        yline(mean(Delta_U_minMax .* delta_half_minMax), 'color', 'green', 'HandleVisibility', 'off')
        yline(mean(Delta_U_middle .* delta_half_middle), 'color', 'blue', 'HandleVisibility', 'off')
    end

end

% Legend
plot(nan, nan, 'linewidth', lw, 'color', 'red', 'displayname', 'Mean')
plot(nan, nan, 'linewidth', lw, 'color', 'green', 'displayname', 'Max/Min')
plot(nan, nan, 'linewidth', lw, 'color', 'blue', 'displayname', 'Centerline')

hold off
xlabel('$x / D$', 'interpreter', 'latex', 'fontsize', 18)
legend('box', 'off', 'location', 'southeast')


%% Show what the different halfwidths look like also


clc; close all
figure('color', 'white')
title('$\delta_{1/2}$', 'interpreter', 'latex', 'fontsize', 24)
hold on
for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(1);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;
    massaged_wake_uncropped = cleaned(block).u;

    % Crop inner side
    raw_wake(Y > 3) = nan;
    massaged_wake(Y > 3) = nan;
    massaged_wake_uncropped(Y > 3) = nan;


    % Plot
    contourf(X, Y, massaged_wake_uncropped, 'linestyle', 'none', 'HandleVisibility', 'off')
    plot(x, halfwidth_mean.massaged(block, :), 'color', 'red', 'linewidth', lw, 'HandleVisibility', 'off')
    plot(x, halfwidth_minMax.massaged(block, :), 'color', 'green', 'linewidth', lw, 'HandleVisibility', 'off')
    plot(x, halfwidth_middle.massaged(block, :), 'color', 'blue', 'linewidth', lw, 'HandleVisibility', 'off')
    
end

% Legend
plot(nan, nan, 'linewidth', lw, 'color', 'red', 'displayname', 'Mean')
plot(nan, nan, 'linewidth', lw, 'color', 'green', 'displayname', 'Max/Min')
plot(nan, nan, 'linewidth', lw, 'color', 'blue', 'displayname', 'Centerline')


hold off
colorbar()
axis equal
set(gca, 'YDir', 'reverse')
ylim([-3, 6])
yline(0, 'linestyle', '--', 'HandleVisibility', 'off')
xlabel('$x / D$', 'interpreter', 'latex')
ylabel('$z / D$', 'interpreter', 'latex')
legend('box', 'off', 'location', 'eastoutside')

%% Test momentum thickness integrand profiles

turbine = 1;
block = 3;

Y = data(block).Y - turbine_positions(turbine);
y = Y(:,1);

tmp_data = cleaned(block).u;
% tmp_centerline = mean_centerline_velocity(block).massaged.u;

idx = 100;

% Plot profile of data
clc; close all
figure('color', 'white')
hold on


% Momentum thickness integrand
u = tmp_data(:, idx);
% u1 = max(u, [], 'all', 'omitnan');
u1 = outer_mean(block).u(idx);

% u2 = min(u, [], 'all', 'omitnan');
u2 = inner_mean(block).u(idx);

mom_integrand = ((u - u2) .* (u1 - u)) / ((u1 - u2)^2);
plot(mom_integrand, y)

hold off
xline(0)
set(gca, 'YDir', 'reverse')
xlim([-1, 1])

clear Y y tmp_data ids u u1 u2 mom_integrand turbine block u_inf

%% Loop through to better check

turbine = 1;

% Near, middle, and far regions
colors(1).c = 'red';
colors(2).c = 'green';
colors(3).c = 'blue';

skip = 10;

figure('color', 'white')
hold on
for block = 2:3

    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    tmp_data = cleaned(block).u;

    for i = 1:skip:size(tmp_data, 2)
        u = tmp_data(:, i);
        % u1 = max(u, [], 'all', 'omitnan');
        % u2 = min(u, [], 'all', 'omitnan');

        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);

        mom_integrand = ((u - u2) .* (u1 - u)) / ((u1 - u2)^2);

        % Hard mask near the mixing layer
        mom_integrand(y > 4.5) = nan;

        % 1. Find the main mixing-layer peak (should be near y ~ 0)
        [~, peak_idx] = max(mom_integrand);
        
        % 2. Search downward from the peak for the first point that drops
        %    below a small threshold (i.e. the integrand has "closed out")
        thresh = 0.001;   % tune: ~2% of the peak value works well
        search = mom_integrand(peak_idx:end);
        below = find(search < thresh, 1, 'first');
        
        if isempty(below)
            cutoff_idx = length(y);   % never closes out within FOV
            disp('xxx')
        else
            cutoff_idx = peak_idx + below - 1;
        end
        
        % 3. Mask everything beyond cutoff_idx
        mom_integrand(cutoff_idx+1:end) = NaN;

        % Plot
        plot(mom_integrand, y, 'color', colors(block).c)
    end

end
hold off
hold off
xline(0)
set(gca, 'YDir', 'reverse')
xlim([-1, 1])
yline(0)
yline(3)

clear below block cutoff_idx i idx mom_integrand peak_idx search skip thresh
clear u u1 u2 y Y turbine


%% Claude test (compute momentum thickness) + MEAN values

turbine = 1;

% Near, middle, and far regions
colors(1).c = 'red';
colors(2).c = 'green';
colors(3).c = 'blue';
skip = 1;
thresh_frac = 0.001;   % fraction of peak below which we consider integrand "closed"

theta = cell(1,3);
x_theta = cell(1,3);

figure('color', 'white')
hold on
for block = 1:3
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);
    X = data(block).X;
    x = X(1,:);
    tmp_data = cleaned(block).u;

    theta{block}   = nan(1, size(tmp_data, 2));
    x_theta{block} = x;

    for i = 1:size(tmp_data, 2)
        u  = tmp_data(:, i);
        % u1 = max(u, [], 'omitnan');
        % u2 = min(u, [], 'omitnan');

        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);
        if ~isfinite(u1) || ~isfinite(u2) || (u1 - u2) < 1e-6
            continue
        end

        mom_integrand = ((u - u2) .* (u1 - u)) / ((u1 - u2)^2);

        % Find main peak then first drop below threshold after it
        [pk, peak_idx] = max(mom_integrand);
        thresh = thresh_frac * pk;
        below  = find(mom_integrand(peak_idx:end) < thresh, 1, 'first');
        if isempty(below)
            cutoff_idx = numel(y);
        else
            cutoff_idx = peak_idx + below - 1;
        end

        % Mask everything beyond the mixing-layer region
        mom_integrand(cutoff_idx+1:end) = NaN;

        % Integrate (only over finite points)
        good = isfinite(mom_integrand) & isfinite(y);
        if nnz(good) >= 3
            theta{block}(i) = trapz(y(good), mom_integrand(good));
        end

        % Plot only the columns we'd otherwise have plotted
        if mod(i-1, skip) == 0
            plot(mom_integrand, y, 'color', colors(block).c)
        end
    end
end
hold off
xline(0)
set(gca, 'YDir', 'reverse')
xlim([-1, 1])
yline(0)
yline(3)

% Quick look at theta(x)
figure('color','white'); hold on
for block = 1:3
    plot(x_theta{block}, theta{block}, 'color', colors(block).c, 'linewidth', 1.5)
end
axis equal
xlabel('x/D')
ylabel('\theta [D]') 
grid on

% Fit a line through blocks 2 and 3
xs = [x_theta{1}, x_theta{2}, x_theta{3}];
ys = [theta{1}, theta{2}, theta{3}];

x_cutoff = 6;
[~, x_cutoff_idx] = min(abs(xs - x_cutoff));
xs_fit = xs(x_cutoff_idx:end);
ys_fit = ys(x_cutoff_idx:end);

P = polyfit(xs_fit, ys_fit, 1);
x_fit = x_cutoff:50;
plot(x_fit, polyval(P, x_fit), 'color', 'black', 'linestyle', '--')
xline(x_cutoff)


clear below block cutoff_idx i idx mom_integrand peak_idx search skip thresh
clear u u1 u2 y Y turbine pk P thresh_frac tmp_data x x)cutoff c_cutoff_idx_x_fit 
clear xs xs_fit X ys


%% Claude test (compute momentum thickness) + MIN/MAX values

turbine = 1;

% Near, middle, and far regions
colors(1).c = 'red';
colors(2).c = 'green';
colors(3).c = 'blue';
skip = 1;
thresh_frac = 0.1;   % fraction of peak below which we consider integrand "closed"

theta = cell(1,3);
x_theta = cell(1,3);

figure('color', 'white')
hold on
for block = 1:3
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);
    X = data(block).X;
    x = X(1,:);
    tmp_data = cleaned(block).u;

    theta{block}   = nan(1, size(tmp_data, 2));
    x_theta{block} = x;

    for i = 1:size(tmp_data, 2)
        u  = tmp_data(:, i);
        u1 = max(u, [], 'omitnan');
        u2 = min(u, [], 'omitnan');

        % u1 = outer_mean(block).u(i);
        % u2 = inner_mean(block).u(i);
        if ~isfinite(u1) || ~isfinite(u2) || (u1 - u2) < 1e-6
            continue
        end

        mom_integrand = ((u - u2) .* (u1 - u)) / ((u1 - u2)^2);

        % Find main peak then first drop below threshold after it
        [pk, peak_idx] = max(mom_integrand);
        thresh = thresh_frac * pk;
        below  = find(mom_integrand(peak_idx:end) < thresh, 1, 'first');
        if isempty(below)
            cutoff_idx = numel(y);
        else
            cutoff_idx = peak_idx + below - 1;
        end

        % Mask everything beyond the mixing-layer region
        mom_integrand(cutoff_idx+1:end) = NaN;

        % Integrate (only over finite points)
        good = isfinite(mom_integrand) & isfinite(y);
        if nnz(good) >= 3
            theta{block}(i) = trapz(y(good), mom_integrand(good));
        end

        % Plot only the columns we'd otherwise have plotted
        if mod(i-1, skip) == 0
            plot(mom_integrand, y, 'color', colors(block).c)
        end
    end
end
hold off
xline(0)
set(gca, 'YDir', 'reverse')
xlim([-1, 1])
yline(0)
yline(3)

% Quick look at theta(x)
figure('color','white'); hold on
for block = 1:3
    plot(x_theta{block}, theta{block}, 'color', colors(block).c, 'linewidth', 1.5)
end
axis equal
xlabel('x/D')
ylabel('\theta [D]') 
grid on

% Fit a line through blocks 2 and 3
xs = [x_theta{1}, x_theta{2}, x_theta{3}];
ys = [theta{1}, theta{2}, theta{3}];

x_cutoff = 6;
[~, x_cutoff_idx] = min(abs(xs - x_cutoff));
xs_fit = xs(x_cutoff_idx:end);
ys_fit = ys(x_cutoff_idx:end);

P = polyfit(xs_fit, ys_fit, 1);
x_fit = x_cutoff:50;
plot(x_fit, polyval(P, x_fit), 'color', 'black', 'linestyle', '--')
xline(x_cutoff)


clear below block cutoff_idx i idx mom_integrand peak_idx search skip thresh
clear u u1 u2 y Y turbine pk P thresh_frac tmp_data x x)cutoff c_cutoff_idx_x_fit 
clear xs xs_fit X ys





%% Try scaling profiles with all the computed values (all blocks overlayed)

turbine = 1;
skip = 20;


figure('color', 'white')
hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth_mean.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195
        u_profile = u(:, i);
        half_width_local = half_width(i);
        momentum_thickness_local = momentum_thickness(i);

        % Local inner/outer flows
        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);

        % Scaled coordinate
        eta_mom = (y - half_width_local) / momentum_thickness_local;

        % Scaled velocity
        u_tilde = (u_profile - u2) / (u1 - u2);

        plot(u_tilde, eta_mom, 'color', colors(block).c)
    end
end

x_plot = -6:0.01:6;
plot(0.5*(1 - erf(sqrt(pi)*x_plot/2)), x_plot, 'color', 'black', 'linewidth', 2)
% plot(0.5*(1 - tanh(x_plot)), x_plot, 'color', 'black')

hold off
axis square
box on
set(gca, 'YDir', 'reverse')

labelFontSize = 18;
xlabel('$\tilde{u}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\eta_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([-4, 4])
% xlim([0, 1])


%% Seperated by tile
figure('color', 'white')
tile = tiledlayout(1,3);

for block = 1:3

    h(block) = nexttile;
    hold on
    title(sprintf('Block %1.0f', block))

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth_mean.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195
        u_profile = u(:, i);
        half_width_local = half_width(i);
        momentum_thickness_local = momentum_thickness(i);

        % Local inner/outer flows
        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);

        % Scaled coordinate
        eta_mom = (y - half_width_local) / momentum_thickness_local;

        % Scaled velocity
        u_tilde = (u_profile - u2) / (u1 - u2);

        plot(u_tilde, eta_mom, 'color', colors(block).c)
    end
    axis square
    set(h(block), 'YDir', 'reverse')

    x_plot = -6:0.01:6;
    plot(0.5*(1 - erf(sqrt(pi)*x_plot/2)), x_plot, 'color', 'black', 'linewidth', 2)

    hold off
end


linkaxes(h, 'xy')

labelFontSize = 18;
xlabel(tile, '$\tilde{u}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$\eta_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylim([-4, 4])
% xlim([0, 1])


%% Self-similarity test for farm-wake edge mixing layer (hub-height, spanwise y)
% Assumes:
% 
%   X = data(block).X;  Y = data(block).Y - turbine_positions(turbine);
%   U = cleaned(block).u  (Ny x Nx)
%   y is in D (rotor diameters), U in m/s
% Outputs:
%   (1) full spanwise profiles U(y) at several x stations
%   (2) collapse plots: U~ vs eta for left and right edges

u_inf  = 8;
turbine = 1;

blocks_to_use = 2:3;

% ---- user-tunable knobs ----
x_stations_per_block = 9;           % how many downstream stations to sample per block
core_band_D = 1;                    % |y| <= core_band_D used to estimate U2(x)
smooth_y_window = 1;                % smoothing in y (odd integer). set 1 to disable
edge_side = "left";                 % "left", "right", or "both"
use_massaged = true;                % true: cleaned(block).u, false: data(block).u

% optional: limit similarity analysis to far wake range in x (in same units as x array)
% set empty [] to use full block range
x_farwake_min = [];   % e.g. 20;   % (in D if x is already in D)
x_farwake_max = [];   % e.g. 50;

% 5–95% cut for defining edge region (helps avoid overshoots/noise dominating)
use_level_window = true;
level_lo = 0.05;
level_hi = 0.95;

% Collect similarity curves
curves_left  = struct('eta',[],'Uhat',[],'x',[],'block',[]);
curves_right = struct('eta',[],'Uhat',[],'x',[],'block',[]);

% Figure 1: full spanwise profiles
% figure('color','white'); tiledlayout(1,3,'padding','compact','tilespacing','compact');
for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U = cleaned(block).u;
    else
        U = data(block).u;
    end

    % inner-side crop (if you want to keep it consistent with your earlier code)
    U(Y > 3) = nan;

    % choose x stations (optionally within far-wake range)
    x_mask = true(size(x));
    if ~isempty(x_farwake_min), x_mask = x_mask & (x >= x_farwake_min); end
    if ~isempty(x_farwake_max), x_mask = x_mask & (x <= x_farwake_max); end
    x_idx_all = find(x_mask);

    if numel(x_idx_all) < x_stations_per_block
        x_idx = x_idx_all;
    else
        x_idx = round(linspace(x_idx_all(1), x_idx_all(end), x_stations_per_block));
    end

    % h(bi) = nexttile(bi); hold on
    % for k = 1:numel(x_idx)
    %     j = x_idx(k);
    %     plot(y, U(:,j), 'linewidth', 1.2, 'DisplayName', sprintf('x=%.1f', x(j)));
    % end
    % yline(u_inf,'k-','u_\infty','Interpreter','tex');
    % xlabel('y/D'); ylabel('U (m/s)');
    % title(sprintf('Block %d: full spanwise profiles', block));
    % grid on; box on
end

% linkaxes(h, 'xy')




% Figure 2: collapse plots (left/right edges) USING mean_centerline_velocity as U2
figure('color','white'); tiledlayout(1,1,'padding','compact','tilespacing','compact');

for bi = 2:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    % keep same crop as you've been using
    U(Y > 3) = nan;

    % choose x stations (optionally within far-wake range)
    x_mask = true(size(x));
    if ~isempty(x_farwake_min), x_mask = x_mask & (x >= x_farwake_min); end
    if ~isempty(x_farwake_max), x_mask = x_mask & (x <= x_farwake_max); end
    x_idx_all = find(x_mask);
    if isempty(x_idx_all), x_idx_all = 1:numel(x); end

    if numel(x_idx_all) < x_stations_per_block
        x_idx = x_idx_all;
    else
        x_idx = round(linspace(x_idx_all(1), x_idx_all(end), x_stations_per_block));
    end

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        % optional smoothing in y
        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- U2(x) FROM mean_centerline_velocity (interpolated onto x(j))
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        U1 = u_inf;
        dU = U1 - U2;
        if ~isfinite(dU) || dU <= 1e-6, continue; end

        % Normalized velocity
        Uhat_full = (Uuse - U2) ./ dU;

        % Choose edges
        sides = {};
        if edge_side == "left"  || edge_side == "both",  sides{end+1} = "left";  end
        if edge_side == "right" || edge_side == "both",  sides{end+1} = "right"; end

        for s = 1:numel(sides)
            side = sides{s};

            if side == "left"
                edge_mask = good & (y <= 0);
            else
                edge_mask = good & (y >= 0);
            end

            % restrict to 5–95% region (recommended)
            if use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end

            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            % sort by y
            [y_edge, sortIdx] = sort(y_edge);
            Uhat = Uhat(sortIdx);

            % y0 at 50%
            [~, idx50] = min(abs(Uhat - 0.5));
            y0 = y_edge(idx50);

            % vorticity thickness: delta_omega = 1 / max|d(Uhat)/dy|
            dUhat_dy = gradient(Uhat, y_edge);
            shear_max = max(abs(dUhat_dy));
            if ~isfinite(shear_max) || shear_max <= 1e-9, continue; end
            delta_omega = 1 / shear_max;

            eta = (y_edge - y0) ./ delta_omega;

            nexttile(1); hold on
            plot(Uhat, eta, 'linewidth', 1.2);
            title('Left edge collapse (vorticity thickness)');
            grid on; box on; axis square
            title('Vorticity Thickness Mixing Layer Normalization');
            ylabel('$\eta_{1/2} = \frac{z - z_{1/2}}{\delta_{\omega}}$', 'interpreter', 'latex');
            xlabel('$\tilde{U} = \frac{u - u_{c}}{u_{\infty} - u_{c}}$', 'interpreter', 'latex');
            set(gca, 'YDir', 'reverse')

        end
    end
end

xlim([0, 1])
ylim([-2, 1])


%% Half-width collapse (USING mean_centerline_velocity as U2)

u_inf  = 8;
turbine = 1;
blocks_to_use = 2:3;

% ---- user-tunable knobs ----
x_stations_per_block = 9;           % how many downstream stations to sample per block
core_band_D = 1;                    % |y| <= core_band_D used to estimate U2(x)
smooth_y_window = 1;                % smoothing in y (odd integer). set 1 to disable
use_massaged = true;                % true: cleaned(block).u, false: data(block).u



figure('color','white'); tiledlayout(1,1,'padding','compact','tilespacing','compact');

for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    U(Y > 3) = nan;

    % choose x stations
    x_mask = true(size(x));
    if ~isempty(x_farwake_min), x_mask = x_mask & (x >= x_farwake_min); end
    if ~isempty(x_farwake_max), x_mask = x_mask & (x <= x_farwake_max); end
    x_idx_all = find(x_mask);
    if isempty(x_idx_all), x_idx_all = 1:numel(x); end

    if numel(x_idx_all) < x_stations_per_block
        x_idx = x_idx_all;
    else
        x_idx = round(linspace(x_idx_all(1), x_idx_all(end), x_stations_per_block));
    end

    % Colors based on x location
    colors = slanCM('viridis', length(x_idx));

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- U2(x) from mean_centerline_velocity
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        U1 = u_inf;
        dU = U1 - U2;
        if ~isfinite(dU) || dU <= 1e-6, continue; end

        Uhat_full = (Uuse - U2) ./ dU;

        % --- core center y_c: where U is minimum (optionally within a band)
        % If you want a band: center_search = good & abs(y)<=1.0; else just use good
        center_search = good;  % or: good & (abs(y) <= 1.0);
        [~, imin] = min(Uuse(center_search));
        y_candidates = y(center_search);
        y_c = y_candidates(imin);

        % Choose edges relative to y_c
        sides = {};
        if edge_side == "left"  || edge_side == "both",  sides{end+1} = "left";  end
        if edge_side == "right" || edge_side == "both",  sides{end+1} = "right"; end

        for s = 1:numel(sides)
            side = sides{s};

            if side == "left"
                edge_mask = good & (y <= y_c);
            else
                edge_mask = good & (y >= y_c);
            end

            if use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end

            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            [y_edge, sortIdx] = sort(y_edge);
            Uhat = Uhat(sortIdx);

            % y_{1/2}: location where Uhat ~ 0.5
            [~, idx50] = min(abs(Uhat - 0.5));
            y_half = y_edge(idx50);

            % half-width thickness
            delta_half = abs(y_half - y_c);
            if ~isfinite(delta_half) || delta_half <= 1e-6, continue; end

            eta_hw = (y_edge - y_half) ./ delta_half;


            label = sprintf('$x / D = %2.1f$', x(j));
            nexttile(1); hold on
            plot(Uhat, eta_hw, 'linewidth', 1.2, 'color', colors(k,:), ...
                 'displayname', label);
            title('Half-width Mixing Layer Normalization');
            ylabel('$\eta_{1/2} = \frac{z - z_{1/2}}{\delta_{1/2}}$', 'interpreter', 'latex');
            xlabel('$\tilde{U} = \frac{u - u_{c}}{u_{\infty} - u_{c}}$', 'interpreter', 'latex');
            grid on; box on; axis square
            ylim([-0.65, 0.65])
            xlim([0, 1])
            set(gca, 'YDir', 'reverse')

        end
    end
    % Split up legend based on blocks
    plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
end

leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', 10);
leg.Layout.Tile = 'east';




%% Momentum-thickness collapse (USING mean_centerline_velocity as U2)

u_inf  = 8;
turbine = 1;
blocks_to_use = 2:3;

x_stations_per_block = 9;
smooth_y_window = 1;
use_massaged = true;

figure('color','white');
tiledlayout(1,1,'padding','compact','tilespacing','compact');

for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U    = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U    = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    % keep your crop
    U(Y > 3) = nan;

    % choose x stations
    x_mask = true(size(x));
    if exist('x_farwake_min','var') && ~isempty(x_farwake_min), x_mask = x_mask & (x >= x_farwake_min); end
    if exist('x_farwake_max','var') && ~isempty(x_farwake_max), x_mask = x_mask & (x <= x_farwake_max); end
    x_idx_all = find(x_mask);
    if isempty(x_idx_all), x_idx_all = 1:numel(x); end

    if numel(x_idx_all) < x_stations_per_block
        x_idx = x_idx_all;
    else
        x_idx = round(linspace(x_idx_all(1), x_idx_all(end), x_stations_per_block));
    end

    colors = slanCM('viridis', length(x_idx));

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- U2(x) from mean_centerline_velocity
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        U1 = u_inf;
        dU = U1 - U2;
        if ~isfinite(dU) || dU <= 1e-6, continue; end

        % normalized velocity
        Uhat_full = (Uuse - U2) ./ dU;

        % --- core center y_c (your existing method)
        center_search = good;
        [~, imin] = min(Uuse(center_search));
        y_candidates = y(center_search);
        y_c = y_candidates(imin);

        % Choose edges relative to y_c
        sides = {};
        if edge_side == "left"  || edge_side == "both",  sides{end+1} = "left";  end
        if edge_side == "right" || edge_side == "both",  sides{end+1} = "right"; end

        for s = 1:numel(sides)
            side = sides{s};

            if side == "left"
                edge_mask = good & (y <= y_c);
            else
                edge_mask = good & (y >= y_c);
            end

            % IMPORTANT for theta: only integrate where 0<=Uhat<=1
            % (this enforces two-stream mixing layer region)
            edge_mask = edge_mask & (Uhat_full >= 0) & (Uhat_full <= 1);

            % optional: also apply your 5–95% window if you want
            if exist('use_level_window','var') && use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end

            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            % sort for integration + locating y0
            [y_edge, sortIdx] = sort(y_edge);
            Uhat = Uhat(sortIdx);

            % y0: location where Uhat ~ 0.5 (same anchor as half-width)
            [~, idx50] = min(abs(Uhat - 0.5));
            y0 = y_edge(idx50);

            % Momentum thickness theta = ∫ Uhat(1-Uhat) dy
            integrand = Uhat .* (1 - Uhat);
            theta = trapz(y_edge, integrand);

            if ~isfinite(theta) || theta <= 1e-8, continue; end

            % similarity coordinate based on theta
            eta_theta = (y_edge - y0) ./ theta;

            label = sprintf('$x / D = %2.1f$', x(j));
            nexttile(1); hold on
            plot(Uhat, eta_theta, 'linewidth', 1.2, 'color', colors(k,:), ...
                 'displayname', label);

            title('Momentum-thickness Mixing Layer Normalization');
            ylabel('$\eta_{\theta} = \frac{y - y_{0}}{\theta}$', 'interpreter', 'latex');
            xlabel('$\tilde{U} = \frac{u - u_{c}}{u_{\infty} - u_{c}}$', 'interpreter', 'latex');
            grid on; box on; axis square

            % axis limits: theta-based eta is often wider than half-width
            ylim([-4, 4])
            xlim([0, 1])
            set(gca, 'YDir', 'reverse')
        end
    end

    % Split up legend based on blocks
    plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
end

leg = legend('interpreter', 'latex', 'box', 'off', 'fontsize', 10);
leg.Layout.Tile = 'east';






%% ============================================================
%  LOCAL OUTER STREAM NORMALIZATION (compare vs constant u_inf)
%  U1(x) estimated from an outer band on the measured edge
%  U2(x) taken from mean_centerline_velocity as before
% ============================================================

u_inf  = 8;                 % keep for reference plots
turbine = 1;
blocks_to_use = 2:3;

x_stations_per_block = 9;
smooth_y_window = 1;
use_massaged = true;

% ---- Choose which measured edge you have ----
edge_side = "left";         % "left" or "right"

% ---- Outer-band definition (tune these) ----
outer_frac = 0.15;          % use outermost 15% of available y-range on that edge
outer_buffer = 0.25;        % exclude points within 0.25D of y_c to avoid wake core influence

% ---- For stability: restrict to the "mixing region" when defining y0/theta ----
use_level_window = true;
level_lo = 0.05;
level_hi = 0.95;

% ---- Figure layout: 2x2 (half-width & theta; constant vs local outer) ----
figure('color','white', 'units', 'centimeters', 'position', [2, 2, 20,20]);
tiledlayout(2,2,'padding','tight','tilespacing','tight');

ax_hw_const  = nexttile(1); hold on; title('Half-width: U_1 = u_\infty');
ax_hw_local  = nexttile(2); hold on; title('Half-width: U_1 = U_{out}(x)');
ax_th_const  = nexttile(3); hold on; title('Momentum thickness: U_1 = u_\infty');
ax_th_local  = nexttile(4); hold on; title('Momentum thickness: U_1 = U_{out}(x)');

for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U    = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U    = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    % match your crop
    U(Y > 3) = nan;

    % choose x stations
    x_idx = round(linspace(1, numel(x), x_stations_per_block));


    if block == 2
        colors = slanCM('viridis', length(x_idx));
    elseif block == 3
        colors = slanCM('magma', length(x_idx));
    end

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- core velocity U2(x) from mean_centerline_velocity
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        % --- core center y_c from minimum U (your existing approach)
        center_search = good;
        [~, imin] = min(Uuse(center_search));
        y_candidates = y(center_search);
        y_c = y_candidates(imin);

        % --- define which edge we keep (you said you only have one)
        if edge_side == "left"
            edge_mask_full = good & (y <= y_c);
        else
            edge_mask_full = good & (y >= y_c);
        end
        if nnz(edge_mask_full) < 8, continue; end

        % --- LOCAL OUTER STREAM Uout(x): outermost band on that edge
        y_edge_all = y(edge_mask_full);
        U_edge_all = Uuse(edge_mask_full);

        % Exclude buffer near y_c
        outer_keep = abs(y_edge_all - y_c) >= outer_buffer;

        % If too few points after buffer, relax buffer
        if nnz(outer_keep) < 5
            outer_keep = true(size(outer_keep));
        end

        y_edge_buff = y_edge_all(outer_keep);
        U_edge_buff = U_edge_all(outer_keep);

        % Outer band: farthest |y - y_c| values on that edge
        dist = abs(y_edge_buff - y_c);
        [dist_sorted, idx_sorted] = sort(dist, 'descend');

        n_outer = max(5, round(outer_frac * numel(dist_sorted)));
        idx_outer = idx_sorted(1:n_outer);

        Uout = median(U_edge_buff(idx_outer), 'omitnan');

        % --- Now build two normalizations:
        % (A) constant outer: U1 = u_inf
        % (B) local outer:    U1 = Uout

        % Helper function: compute half-width eta and theta-eta for a given U1
        % (implemented inline for clarity)

        U1_list = [u_inf, Uout];
        for u1_case = 1:2
            U1 = U1_list(u1_case);
            dU = U1 - U2;
            if ~isfinite(dU) || dU <= 1e-6
                continue
            end

            Uhat_full = (Uuse - U2) ./ dU;

            % Edge mask restricted to 0<=Uhat<=1 for theta stability
            edge_mask = edge_mask_full & (Uhat_full >= 0) & (Uhat_full <= 1);

            if use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end
            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            [y_edge, idxS] = sort(y_edge);
            Uhat = Uhat(idxS);

            % y0 at 50% location
            [~, idx50] = min(abs(Uhat - 0.5));
            y0 = y_edge(idx50);

            % ---- Half-width coordinate (need y_{1/2} and delta_{1/2})
            y_half = y0;                  % same anchor
            delta_half = abs(y_half - y_c);
            if ~isfinite(delta_half) || delta_half <= 1e-6
                continue
            end
            eta_hw = (y_edge - y_half) ./ delta_half;

            % ---- Momentum thickness theta and eta_theta
            theta = trapz(y_edge, Uhat .* (1 - Uhat));
            if ~isfinite(theta) || theta <= 1e-8
                continue
            end
            eta_th = (y_edge - y0) ./ theta;

            % ---- Plot (match your style: Uhat on x-axis, eta on y-axis)
            label = sprintf('$x/D = %2.1f$', x(j));

            if u1_case == 1
                plot(ax_hw_const, Uhat, eta_hw, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
                plot(ax_th_const, Uhat, eta_th, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
            else
                plot(ax_hw_local, Uhat, eta_hw, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
                plot(ax_th_local, Uhat, eta_th, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
            end
        end
    end

    % spacer for legend grouping
    plot(ax_hw_const, nan, nan, 'color','white','DisplayName',' ');
    plot(ax_hw_local, nan, nan, 'color','white','DisplayName',' ');
    plot(ax_th_const, nan, nan, 'color','white','DisplayName',' ');
    plot(ax_th_local, nan, nan, 'color','white','DisplayName',' ');
end

% Cosmetics
axes_list = [ax_hw_const, ax_hw_local, ax_th_const, ax_th_local];
for ax = axes_list
    grid(ax,'on'); box(ax,'on'); axis(ax,'square');
    set(ax, 'YDir', 'reverse');
    xlim(ax, [0 1]);
end

ylabel(ax_hw_const, '$\eta_{1/2}$', 'interpreter','latex');
ylabel(ax_hw_local, '$\eta_{1/2}$', 'interpreter','latex');
ylabel(ax_th_const, '$\eta_{\theta}$', 'interpreter','latex');
ylabel(ax_th_local, '$\eta_{\theta}$', 'interpreter','latex');

xlabel(ax_hw_const, '$\tilde{U}$', 'interpreter','latex');
xlabel(ax_hw_local, '$\tilde{U}$', 'interpreter','latex');
xlabel(ax_th_const, '$\tilde{U}$', 'interpreter','latex');
xlabel(ax_th_local, '$\tilde{U}$', 'interpreter','latex');

% Reasonable y-lims (tune after first run)
ylim(ax_hw_const, [-0.65 0.65]);
ylim(ax_hw_local, [-0.65 0.65]);
ylim(ax_th_const, [-4 4]);
ylim(ax_th_local, [-4 4]);

% One legend (use the local half-width panel; you can move it)
leg = legend(ax_hw_local, 'interpreter','latex', 'box','off', 'fontsize', 10);
leg.Layout.Tile = 'east';

linkaxes([ax_hw_const, ax_hw_local], 'xy')
linkaxes([ax_th_const, ax_th_local], 'xy')









%% ============================================================
%  LOCAL OUTER STREAM NORMALIZATION (compare vs constant u_inf)
%  U1(x) estimated from an outer band on the measured edge
%  U2(x) taken from mean_centerline_velocity as before

% VELOCITY PROFILES
% ============================================================

u_inf  = 8;                 % keep for reference plots
turbine = 1;
blocks_to_use = 2:3;

x_stations_per_block = 10;
smooth_y_window = 1;
use_massaged = false;

% ---- Choose which measured edge you have ----
edge_side = "left";         % "left" or "right"

% ---- Outer-band definition (tune these) ----
outer_frac = 0.15;          % use outermost 15% of available y-range on that edge
outer_buffer = 0.25;        % exclude points within 0.25D of y_c to avoid wake core influence

% ---- For stability: restrict to the "mixing region" when defining y0/theta ----
use_level_window = true;
level_lo = 0.01;
level_hi = 0.99;


% Collect all points for global fit (u1_case==2 only, per your current plotting)
U_all   = [];
eta_all = [];


% ---- Figure layout: 2x2 (half-width & theta; constant vs local outer) ----
clc; close all
figure('color','white', 'units', 'centimeters', 'position', [1, 1, 20,20]);
tiledlayout(1,1,'padding','loose','tilespacing','compact');

ax = nexttile(1); hold(ax,'on');          % collapse plot
% ax_theta = nexttile(2); hold(ax_theta,'on');
% ax_hw    = nexttile(3); hold(ax_hw,'on');

% storage (local-outer only)
theta_x_all = [];  theta_all = [];
hw_x_all    = [];  hw_all    = [];
yc_all      = [];  y0_all    = [];  % optional sanity
block_all   = [];


hold on
for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U    = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U    = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    % match your crop
    U(Y > 3) = nan;

    % choose x stations
    x_idx = round(linspace(1, numel(x), x_stations_per_block));


    if block == 2
        colors = slanCM('Reds', length(x_idx));
    elseif block == 3
        colors = slanCM('Blues', length(x_idx));
    end

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- core velocity U2(x) from mean_centerline_velocity
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        % --- core center y_c from minimum U (your existing approach)
        center_search = good;
        [~, imin] = min(Uuse(center_search));
        y_candidates = y(center_search);
        y_c = y_candidates(imin);

        % --- define which edge we keep (you said you only have one)
        if edge_side == "left"
            edge_mask_full = good & (y <= y_c);
        else
            edge_mask_full = good & (y >= y_c);
        end
        if nnz(edge_mask_full) < 8, continue; end

        % --- LOCAL OUTER STREAM Uout(x): outermost band on that edge
        y_edge_all = y(edge_mask_full);
        U_edge_all = Uuse(edge_mask_full);

        % Exclude buffer near y_c
        outer_keep = abs(y_edge_all - y_c) >= outer_buffer;

        % If too few points after buffer, relax buffer
        if nnz(outer_keep) < 5
            outer_keep = true(size(outer_keep));
        end

        y_edge_buff = y_edge_all(outer_keep);
        U_edge_buff = U_edge_all(outer_keep);

        % Outer band: farthest |y - y_c| values on that edge
        dist = abs(y_edge_buff - y_c);
        [dist_sorted, idx_sorted] = sort(dist, 'descend');

        n_outer = max(5, round(outer_frac * numel(dist_sorted)));
        idx_outer = idx_sorted(1:n_outer);

        Uout = median(U_edge_buff(idx_outer), 'omitnan');

        % --- Now build two normalizations:
        % (A) constant outer: U1 = u_inf
        % (B) local outer:    U1 = Uout

        % Helper function: compute half-width eta and theta-eta for a given U1
        % (implemented inline for clarity)

        U1_list = [u_inf, Uout];
        for u1_case = 1:2
            U1 = U1_list(u1_case);
            dU = U1 - U2;
            if ~isfinite(dU) || dU <= 1e-6
                continue
            end

            Uhat_full = (Uuse - U2) ./ dU;

            % Edge mask restricted to 0<=Uhat<=1 for theta stability
            edge_mask = edge_mask_full & (Uhat_full >= 0) & (Uhat_full <= 1);

            if use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end
            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            [y_edge, idxS] = sort(y_edge);
            Uhat = Uhat(idxS);

            % y0 at 50% location
            [~, idx50] = min(abs(Uhat - 0.5));
            y0 = y_edge(idx50);

            % ---- Half-width coordinate (need y_{1/2} and delta_{1/2})
            y_half = y0;                  % same anchor
            delta_half = abs(y_half - y_c);
            if ~isfinite(delta_half) || delta_half <= 1e-6
                continue
            end
            eta_hw = (y_edge - y_half) ./ delta_half;

            % ---- Momentum thickness theta and eta_theta
            theta = trapz(y_edge, Uhat .* (1 - Uhat));
            if ~isfinite(theta) || theta <= 1e-8
                continue
            end
            eta_th = (y_edge - y0) ./ theta;

            % --- Save thickness metrics for local-outer case only
            if u1_case == 2
                theta_x_all = [theta_x_all; x(j)];
                theta_all   = [theta_all; theta];
            
                hw_x_all    = [hw_x_all; x(j)];
                hw_all      = [hw_all; delta_half];
            
                yc_all      = [yc_all; y_c];
                y0_all      = [y0_all; y0];
            
                block_all   = [block_all; block];
            end



            % ---- Plot (match your style: Uhat on x-axis, eta on y-axis)
            label = sprintf('$x/D = %2.1f$', x(j));

            if u1_case == 2
                sz = 10;
                % plot(ax, Uhat, eta_th, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
                scatter(ax, Uhat, eta_th, sz, 'o', 'filled', 'MarkerFaceColor', colors(k,:), 'displayname', label);

                % Save points for global line fit
                U_all   = [U_all;   Uhat(:)];
                eta_all = [eta_all; eta_th(:)];
            end
        end
    end

    % spacer for legend grouping
    plot(ax, nan, nan, 'color','white','DisplayName',' ');
end


% --- Fit a single line to ALL plotted points ---
good_fit = isfinite(U_all) & isfinite(eta_all);
U_fit = U_all(good_fit);
eta_fit = eta_all(good_fit);



line_color = 'black';
line_width = 3;
line_style = ':';
fit_buff = 0.05;

if numel(U_fit) >= 2
    % p(1)=slope, p(2)=intercept

    % Fit a single line to all data
    p = polyfit(U_fit, eta_fit, 1);     
    U_line = linspace(0, 1, 200);
    eta_line = polyval(p, U_line);

    % P = plot(ax, U_line, eta_line, 'k-', 'linewidth', 2, ...
    %     'DisplayName', sprintf('$\\eta_\\theta = %.2f \\cdot \\tilde{u} %+0.2f$', p(1), p(2)));
    % P.Color(4) = 0.25;


    % Break up linear fit into 3 sections
    % Section 1: 0 - 0.25
    sec_fit_mask = (U_fit < 0.25);

    p = polyfit(U_fit(sec_fit_mask), eta_fit(sec_fit_mask), 1);     
    U_line = linspace(0 + fit_buff, 0.25 - fit_buff, 5);
    eta_line = polyval(p, U_line);

    plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
    plot(ax, U_line, eta_line, 'linestyle', line_style, 'linewidth', line_width, 'color', line_color, ...
        'DisplayName', sprintf('$0 < \\tilde{u} < 0.25$\n$\\eta_\\theta = %.2f \\cdot \\tilde{u} %+0.2f$', p(1), p(2)));


    % Section 1: 0.25 - 0.75
    sec_fit_mask = (U_fit > 0.25) & (U_fit < 0.75);

    p = polyfit(U_fit(sec_fit_mask), eta_fit(sec_fit_mask), 1);     
    U_line = linspace(0.25 + fit_buff, 0.75 - fit_buff, 5);
    eta_line = polyval(p, U_line);

    plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
    plot(ax, U_line, eta_line, 'linestyle', line_style, 'linewidth', line_width, 'color', line_color, ...
        'DisplayName', sprintf('$0.25 < \\tilde{u} < 0.75$\n$\\eta_\\theta = %.2f \\cdot \\tilde{u} %+0.2f$', p(1), p(2)));


    % Section 3: .75 - 1
    sec_fit_mask = (U_fit > 0.75);

    p = polyfit(U_fit(sec_fit_mask), eta_fit(sec_fit_mask), 1);     
    U_line = linspace(0.75 + fit_buff, 1 - fit_buff, 5);
    eta_line = polyval(p, U_line);

    plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
    plot(ax, U_line, eta_line, 'linestyle', line_style, 'linewidth', line_width, 'color', line_color, ...
        'DisplayName', sprintf('$0.75 < \\tilde{u} < 1$\n$\\eta_\\theta = %.2f \\cdot \\tilde{u} %+0.2f$', p(1), p(2)));

    % If you want R^2:
    eta_pred = polyval(p, U_fit);
    SSres = sum((eta_fit - eta_pred).^2);
    SStot = sum((eta_fit - mean(eta_fit)).^2);
    R2 = 1 - SSres/SStot;
    fprintf('Global line fit: slope = %.6f, intercept = %.6f, R^2 = %.4f\n', p(1), p(2), R2);
end

xline([0.25, 0.75], 'HandleVisibility', 'off')



% Mixing layer theoretical (manual)
a = -0.5;
eta_tanh = linspace(-4,4,100);
u_tanh = 0.5 * (1 + tanh(a * eta_tanh));
% plot(ax, u_tanh, eta_tanh, 'color', 'magenta', 'linewidth', 4);

% Fit eta = (1/a)*atanh(2u-1) + b
mask = (U_fit > 0.1) & (U_fit < 0.9);
u_m = U_fit(mask);
e_m = eta_fit(mask);

X = atanh(2*u_m - 1);          % predictor
P = polyfit(X, e_m, 1);        % e ~= P1*X + P2
A = P(1);                      % A = 1/a
b = P(2);
a_fit = 1/A;

fprintf('atanh fit: a = %.4f, b = %.4f\n', a_fit, b);

u_line = linspace(0.02, 0.98, 400);
eta_model = (1/a_fit)*atanh(2*u_line - 1) + b;
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
plot(ax, u_line, eta_model, 'm-', 'linewidth', 3, 'DisplayName', 'atanh fit');


hold off


% % ---- Plot theta(x) and half-width delta_{1/2}(x) (local outer case) ----
% 
% % --- theta plot
% good_th = isfinite(theta_x_all) & isfinite(theta_all);
% x_th = theta_x_all(good_th);
% th   = theta_all(good_th);
% blk  = block_all(good_th);
% 
% [x_th, iS] = sort(x_th);
% th  = th(iS);
% blk = blk(iS);
% 
% if ~isempty(x_th)
%     scatter(ax_theta, x_th, th, 25, 'k', 'filled', 'DisplayName', '$\theta$');
%     % plot(ax_theta, x_th, th, 'k-', 'linewidth', 1.2, 'HandleVisibility','off');
% 
%     % linear fit (optional)
%     pth = polyfit(x_th, th, 1);
%     x_line = linspace(min(x_th), max(x_th), 200);
%     th_line = polyval(pth, x_line);
%     % plot(ax_theta, x_line, th_line, 'r-', 'linewidth', 2, ...
%     %     'DisplayName', sprintf('$\\theta = %.3f\\,x %+0.3f$', pth(1), pth(2)));
% end
% 
% grid(ax_theta,'on'); box(ax_theta,'on');
% ylim(ax_theta, [0,1])
% xlabel(ax_theta, '$x/D$', 'interpreter','latex','fontsize',FS);
% ylabel(ax_theta, '$\theta$ (in $D$)', 'interpreter','latex','fontsize',FS);
% title(ax_theta, 'Momentum thickness growth (local $U_1$)','fontsize',FS);
% legend(ax_theta,'interpreter','latex','box','off','location','best');
% 
% % --- half-width plot
% good_hw = isfinite(hw_x_all) & isfinite(hw_all);
% x_hw = hw_x_all(good_hw);
% hw   = hw_all(good_hw);
% 
% [x_hw, iS2] = sort(x_hw);
% hw = hw(iS2);
% 
% if ~isempty(x_hw)
%     scatter(ax_hw, x_hw, hw, 25, 'b', 'filled', 'DisplayName', '$\delta_{1/2}$');
%     % plot(ax_hw, x_hw, hw, 'b-', 'linewidth', 1.2, 'HandleVisibility','off');
% 
%     % linear fit (optional)
%     phw = polyfit(x_hw, hw, 1);
%     x_line = linspace(min(x_hw), max(x_hw), 200);
%     hw_line = polyval(phw, x_line);
%     % plot(ax_hw, x_line, hw_line, 'r-', 'linewidth', 2, ...
%     %     'DisplayName', sprintf('$\\delta_{1/2} = %.3f\\,x %+0.3f$', phw(1), phw(2)));
% end
% 
% grid(ax_hw,'on'); box(ax_hw,'on');
% % set(ax_hw, 'YDir', 'reverse')
% xlabel(ax_hw, '$x/D$', 'interpreter','latex','fontsize',FS);
% ylabel(ax_hw, '$\delta_{1/2}$ (in $D$)', 'interpreter','latex','fontsize',FS);
% title(ax_hw, 'Half-width growth (local $U_1$)','fontsize',FS);
% legend(ax_hw,'interpreter','latex','box','off','location','best');
% 
% 
% % Cosmetics for theta axes
% grid(ax_theta,'on'); box(ax_theta,'on');
% xlabel(ax_theta, '$x/D$', 'interpreter','latex','fontsize',FS);
% ylabel(ax_theta, '$\theta/D$', 'interpreter','latex','fontsize',FS);
% title(ax_theta, 'Momentum thickness growth (local outer normalization)','fontsize',FS);


% Cosmetics
grid(ax,'on'); box(ax,'on'); axis(ax,'square');
set(ax, 'YDir', 'reverse');
xlim(ax, [0 1]);
FS = 18;
ylabel(ax, '$\eta_{\theta}$', 'interpreter','latex','fontsize',FS);
xlabel(ax, '$\tilde{u}$', 'interpreter','latex','fontsize',FS);


% Reasonable y-lims (tune after first run)
ylim([-0.65 0.65]);
ylim([-4 4]);

% One legend (use the local half-width panel; you can move it)
leg = legend('interpreter','latex', 'box','off', 'fontsize', 10);
leg.Layout.Tile = 'east';








%% ============================================================
% MIXING LAYER TANH FIT
% ============================================================

u_inf  = 8;                 % keep for reference plots
turbine = 1;
blocks_to_use = 2:3;

x_stations_per_block = 10;
smooth_y_window = 1;
use_massaged = false;

% ---- Choose which measured edge you have ----
edge_side = "left";         % "left" or "right"

% ---- Outer-band definition (tune these) ----
outer_frac = 0.15;          % use outermost 15% of available y-range on that edge
outer_buffer = 0.25;        % exclude points within 0.25D of y_c to avoid wake core influence

% ---- For stability: restrict to the "mixing region" when defining y0/theta ----
use_level_window = true;
level_lo = 0.01;
level_hi = 0.99;


% Collect all points for global fit (u1_case==2 only, per your current plotting)
U_all   = [];
eta_all = [];
xpt_all = [];   % x/D associated with each (Uhat, eta_th) point



% ---- Figure layout: 2x2 (half-width & theta; constant vs local outer) ----
clc; close all
figure('color','white', 'units', 'centimeters', 'position', [1, 1, 20,20]);
tiledlayout(1,1,'padding','loose','tilespacing','compact');

ax = nexttile(1); hold(ax,'on');          % collapse plot
% ax_theta = nexttile(2); hold(ax_theta,'on');
% ax_hw    = nexttile(3); hold(ax_hw,'on');

% storage (local-outer only)
theta_x_all = [];  theta_all = [];
hw_x_all    = [];  hw_all    = [];
yc_all      = [];  y0_all    = [];  % optional sanity
block_all   = [];


hold on
for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    if use_massaged
        U    = cleaned(block).u;
        U2_x = mean_centerline_velocity(block).massaged.x;
        U2_u = mean_centerline_velocity(block).massaged.u;
    else
        U    = data(block).u;
        U2_x = mean_centerline_velocity(block).raw.x;
        U2_u = mean_centerline_velocity(block).raw.u;
    end

    % match your crop
    U(Y > 3) = nan;

    % choose x stations
    x_idx = round(linspace(1, numel(x), x_stations_per_block));


    if block == 2
        colors = slanCM('Reds', length(x_idx));
    elseif block == 3
        colors = slanCM('Blues', length(x_idx));
    end

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol = U(:,j);
        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) < 10, continue; end

        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % --- core velocity U2(x) from mean_centerline_velocity
        U2 = interp1(U2_x, U2_u, x(j), 'linear', 'extrap');

        % --- core center y_c from minimum U (your existing approach)
        center_search = good;
        [~, imin] = min(Uuse(center_search));
        y_candidates = y(center_search);
        y_c = y_candidates(imin);

        % --- define which edge we keep (you said you only have one)
        if edge_side == "left"
            edge_mask_full = good & (y <= y_c);
        else
            edge_mask_full = good & (y >= y_c);
        end
        if nnz(edge_mask_full) < 8, continue; end

        % --- LOCAL OUTER STREAM Uout(x): outermost band on that edge
        y_edge_all = y(edge_mask_full);
        U_edge_all = Uuse(edge_mask_full);

        % Exclude buffer near y_c
        outer_keep = abs(y_edge_all - y_c) >= outer_buffer;

        % If too few points after buffer, relax buffer
        if nnz(outer_keep) < 5
            outer_keep = true(size(outer_keep));
        end

        y_edge_buff = y_edge_all(outer_keep);
        U_edge_buff = U_edge_all(outer_keep);

        % Outer band: farthest |y - y_c| values on that edge
        dist = abs(y_edge_buff - y_c);
        [dist_sorted, idx_sorted] = sort(dist, 'descend');

        n_outer = max(5, round(outer_frac * numel(dist_sorted)));
        idx_outer = idx_sorted(1:n_outer);

        Uout = median(U_edge_buff(idx_outer), 'omitnan');

        % --- Now build two normalizations:
        % (A) constant outer: U1 = u_inf
        % (B) local outer:    U1 = Uout

        % Helper function: compute half-width eta and theta-eta for a given U1
        % (implemented inline for clarity)

        U1_list = [u_inf, Uout];
        for u1_case = 1:2
            U1 = U1_list(u1_case);
            dU = U1 - U2;
            if ~isfinite(dU) || dU <= 1e-6
                continue
            end

            Uhat_full = (Uuse - U2) ./ dU;

            % Edge mask restricted to 0<=Uhat<=1 for theta stability
            edge_mask = edge_mask_full & (Uhat_full >= 0) & (Uhat_full <= 1);

            if use_level_window
                edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
            end
            if nnz(edge_mask) < 8, continue; end

            y_edge = y(edge_mask);
            Uhat   = Uhat_full(edge_mask);

            [y_edge, idxS] = sort(y_edge);
            Uhat = Uhat(idxS);

            % y0 at 50% location
            [~, idx50] = min(abs(Uhat - 0.5));
            y0 = y_edge(idx50);

            % ---- Half-width coordinate (need y_{1/2} and delta_{1/2})
            y_half = y0;                  % same anchor
            delta_half = abs(y_half - y_c);
            if ~isfinite(delta_half) || delta_half <= 1e-6
                continue
            end
            eta_hw = (y_edge - y_half) ./ delta_half;

            % ---- Momentum thickness theta and eta_theta
            theta = trapz(y_edge, Uhat .* (1 - Uhat));
            if ~isfinite(theta) || theta <= 1e-8
                continue
            end
            eta_th = (y_edge - y0) ./ theta;

            % --- Save thickness metrics for local-outer case only
            if u1_case == 2
                theta_x_all = [theta_x_all; x(j)];
                theta_all   = [theta_all; theta];
            
                hw_x_all    = [hw_x_all; x(j)];
                hw_all      = [hw_all; delta_half];
            
                yc_all      = [yc_all; y_c];
                y0_all      = [y0_all; y0];
            
                block_all   = [block_all; block];
            end



            % ---- Plot (match your style: Uhat on x-axis, eta on y-axis)
            label = sprintf('$x/D = %2.1f$', x(j));

            if u1_case == 2
                sz = 10;
                % plot(ax, Uhat, eta_th, 'linewidth', 1.2, 'color', colors(k,:), 'displayname', label);
                scatter(ax, Uhat, eta_th, sz, 'o', 'filled', 'MarkerFaceColor', colors(k,:), 'displayname', label);

                % Save points for global line fit
                U_all   = [U_all;   Uhat(:)];
                eta_all = [eta_all; eta_th(:)];
                xpt_all = [xpt_all; repmat(x(j), numel(Uhat), 1)];

            end
        end
    end

    % spacer for legend grouping
    plot(ax, nan, nan, 'color','white','DisplayName',' ');
end


% --- Fit a single line to ALL plotted points ---
good_fit = isfinite(U_all) & isfinite(eta_all) & isfinite(xpt_all);
U_fit   = U_all(good_fit);
eta_fit = eta_all(good_fit);
x_fit   = xpt_all(good_fit);



% --- atanh fits in two streamwise bands ---
% model: eta = (1/a)*atanh(2u-1) + b
u_line = linspace(0.02, 0.98, 400);

bands = [21 30; 41 50];
band_names = {'atanh fit (21-30D)','atanh fit (41-50D)'};

plot(nan, nan, 'color', 'white', 'DisplayName', ' ')

for bb = 1:size(bands,1)
    xlo = bands(bb,1);
    xhi = bands(bb,2);

    % band mask + u-cropping for stability
    mask = (x_fit >= xlo) & (x_fit <= xhi) & (U_fit > 0.1) & (U_fit < 0.9);

    if nnz(mask) < 50
        fprintf('Band %d has too few points (%d). Skipping.\n', bb, nnz(mask));
        continue
    end

    u_m = U_fit(mask);
    e_m = eta_fit(mask);

    X = atanh(2*u_m - 1);
    P = polyfit(X, e_m, 1);   % e ~= A*X + b

    A = P(1);
    b = P(2);
    a_fit = 1/A;

    fprintf('%s: a = %.4f, b = %.4f, N = %d\n', band_names{bb}, a_fit, b, nnz(mask));

    eta_model = (1/a_fit)*atanh(2*u_line - 1) + b;

    plot(ax, u_line, eta_model, '-', 'linewidth', 3, ...
        'DisplayName', sprintf('%s: $a=%.3f,\\ b=%.3f$', band_names{bb}, a_fit, b));
end


hold off


% Cosmetics
grid(ax,'on'); box(ax,'on'); axis(ax,'square');
set(ax, 'YDir', 'reverse');
xlim(ax, [0 1]);
FS = 18;
ylabel(ax, '$\eta_{\theta}$', 'interpreter','latex','fontsize',FS);
xlabel(ax, '$\tilde{u}$', 'interpreter','latex','fontsize',FS);


% Reasonable y-lims (tune after first run)
ylim([-0.65 0.65]);
ylim([-4 4]);

% One legend (use the local half-width panel; you can move it)
leg = legend('interpreter','latex', 'box','off', 'fontsize', 10);
leg.Layout.Tile = 'east';















%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TURBULENCE COLLAPSE ON SINGLE FARM EDGE (half-width coords)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Settings (match your earlier choices)
use_massaged = false;
edge_side = "left";           % you said you only have one edge
x_stations_per_block = 5;    % increase a bit for smoother story
smooth_y_window = 3;          % small smoothing helps derivatives
use_level_window = true;
level_lo = 0.05;
level_hi = 0.95;

u_inf = 8;
turbine = 1;
blocks_to_use = 2:3;   % your far wake blocks

figure('color','white');
tiledlayout(2,2,'padding','compact','tilespacing','compact');

axU  = nexttile; hold on; title('Mean collapse');               xlabel('\eta_{1/2}'); ylabel('$\tilde{U}$', 'interpreter', 'latex');
axUV = nexttile; hold on; title('-uv / \DeltaU^2');             xlabel('\eta_{1/2}'); ylabel('$-\overline{u''v''}/\Delta U^2$', 'interpreter', 'latex');
axK  = nexttile; hold on; title('uu, vv / \DeltaU^2');          xlabel('\eta_{1/2}'); ylabel('$\overline{u''^2},\overline{v''^2} / \Delta U^2$', 'interpreter', 'latex');
axP  = nexttile; hold on; title('Production proxy collapse');   xlabel('\eta_{1/2}'); ylabel('$P / (\Delta U^3/\delta_{1/2})$', 'interpreter', 'latex');


for bi = 1:numel(blocks_to_use)
    block = blocks_to_use(bi);

    % Coordinates
    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    x = X(1,:);
    y = Y(:,1);

    % Fields
    if use_massaged
        U  = cleaned(block).u;
        UU = cleaned(block).uu;
        VV = cleaned(block).vv;
        UV = cleaned(block).uv;

        U2u = mean_centerline_velocity(block).massaged.u;  % 1 x Nx (same x grid)
    else
        U  = data(block).u;
        UU = data(block).uu;
        VV = data(block).vv;
        UV = data(block).uv;

        U2u = mean_centerline_velocity(block).raw.u;
    end

    % Crop inner side like you've been doing (keep consistent)
    y_crop = 6;
    U(Y > y_crop)  = nan;  UU(Y > y_crop) = nan;  VV(Y > y_crop) = nan;  UV(Y > y_crop) = nan;

    % Pick x stations (spread over available range)
    x_idx = round(linspace(1, numel(x), x_stations_per_block));

    for k = 1:numel(x_idx)
        j = x_idx(k);

        Ucol  = U(:,j);
        uucol = filloutliers(UU(:,j), 'previous');
        vvcol = filloutliers(VV(:,j), 'previous');
        uvcol = filloutliers(UV(:,j), 'previous');

        good = isfinite(Ucol) & isfinite(y) & isfinite(uucol) & isfinite(vvcol) & isfinite(uvcol);
        if nnz(good) < 15, continue; end

        % Smooth mean profile in y (helps dU/dy)
        Uuse = Ucol;
        if smooth_y_window > 1
            Uuse(good) = movmean(Uuse(good), smooth_y_window, 'omitnan');
        end

        % Core + outer velocities
        U2 = U2u(j);         % <-- THIS IS YOUR REQUEST: use mean_centerline_velocity as U2
        U1 = u_inf;
        dU = U1 - U2;
        if ~isfinite(dU) || dU <= 1e-6, continue; end

        % Normalized mean
        Uhat_full = (Uuse - U2) ./ dU;

        % Single edge mask (left edge only)
        if edge_side == "left"
            edge_mask = good & (y <= 0);
        else
            edge_mask = good & (y >= 0);
        end

        % Optional: only use the actual shear-layer band (5–95%)
        if use_level_window
            edge_mask = edge_mask & (Uhat_full >= level_lo) & (Uhat_full <= level_hi);
        end
        if nnz(edge_mask) < 10, continue; end

        % Extract + sort
        y_edge = y(edge_mask);
        Uhat   = Uhat_full(edge_mask);
        uu     = uucol(edge_mask);
        vv     = vvcol(edge_mask);
        uv     = uvcol(edge_mask);
        Uedge  = Uuse(edge_mask);

        [y_edge, idx] = sort(y_edge);
        Uhat  = Uhat(idx);
        uu    = uu(idx);
        vv    = vv(idx);
        uv    = uv(idx);
        Uedge = Uedge(idx);

        % Half-width location y_{1/2} from Uhat ~ 0.5
        [~, idx50] = min(abs(Uhat - 0.5));
        y_half = y_edge(idx50);

        % Core center y_c: use y-location of minimum mean within edge-side neighborhood
        % (keeps you robust even if wake drifts a bit)
        [~, imin] = min(Uuse(good));
        y_good = y(good);
        y_c = y_good(imin);

        delta_half = abs(y_half - y_c);
        if ~isfinite(delta_half) || delta_half <= 1e-6, continue; end

        eta = (y_edge - y_half) ./ delta_half;

        % Normalize turbulence
        uvN = -uv ./ (dU.^2);                 % Reynolds shear stress (sign flipped for mixing intuition)
        uuN =  uu ./ (dU.^2);
        vvN =  vv ./ (dU.^2);

        % Production proxy: P ~ -uv * dU/dy
        dUdy = gradient(Uedge,  0.08 * y_edge);
        P = (-uv) .* dUdy;                    % units: (m^2/s^2) * (1/s) = m^2/s^3
        PN = P ./ (dU.^3 ./ delta_half);      % dimensionless-ish production scaling

        % Plot
        plot(axU,  eta, Uhat, 'linewidth', 1.2);
        plot(axUV, eta, uvN,  'linewidth', 1.2);
        plot(axK,  eta, uuN,  'linewidth', 1.2);
        plot(axK,  eta, vvN,  '--',        'linewidth', 1.2);
        plot(axP,  eta, PN,   'linewidth', 1.2);
    end
end

% Cosmetics
grid(axU,'on');  box(axU,'on');  ylim(axU,[0 1]); xlim(axU,[-2 2]);
grid(axUV,'on'); box(axUV,'on'); xlim(axUV,[-2 2]);
grid(axK,'on');  box(axK,'on');  xlim(axK,[-2 2]);
grid(axP,'on');  box(axP,'on');  xlim(axP,[-2 2]);
legend(axK, {'uu/\DeltaU^2','vv/\DeltaU^2'}, 'location','best');
