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
clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);

ax = gca;
ax.FontSize = 8;
set(gca, 'TickLabelInterpreter', 'latex');

hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y, cleaned(block).u / 8, 100, 'linestyle', 'none')
    clim([0.2, 1])
end
hold off
colormap(slanCM('gothic'))
axis equal
set(gca, 'YDir', 'reverse')
% xlabel('$u \mathbin{/} u_{\infty}$', 'interpreter', 'latex')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')

xlim([0, 50])
ylim([-12, 4.5])
yticks(-12:3:3)
% C = colorbar();

% xlabel('$u \mathbin{/} u_{\infty}$', 'interpreter', 'latex')
C = colorbar;
C.Label.String = '$u \mathbin{/} u_{\infty}$';
C.Label.Interpreter = 'latex';
C.Label.FontSize = 10;

% pause(3)
% exportgraphics(fig, fullfile('/Users/zeinsadek/Downloads', 'NAWEA_Farm_tmp.pdf'), 'Resolution', 600)


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
    plot(xs(block).x, outer_mean(block).u, 'color', 'black', 'linewidth', lw, 'HandleVisibility', 'off')

    % Min/max
    plot(xs(block).x, inner_min(block).u, 'color', 'red', 'linestyle', '--', 'linewidth', lw, 'HandleVisibility', 'off')
    plot(xs(block).x, outer_max(block).u, 'color', 'black', 'linestyle', '--', 'linewidth', lw, 'HandleVisibility', 'off')

    % Centerline of middle farm
    plot(xs(block).x, middle_centerline(block).u, 'color', 'blue', 'HandleVisibility', 'off', 'linewidth', lw)

end

% Legend (colors)
plot(nan, nan, 'linewidth', lw, 'color', 'black', 'DisplayName', 'Outer Region')
plot(nan, nan, 'linewidth', lw, 'color', 'red', 'DisplayName', 'Inner Region')
plot(nan, nan, 'linewidth', lw, 'color', 'white', 'DisplayName', ' ')
plot(nan, nan, 'linewidth', lw, 'color', 'black', 'linestyle', '-', 'DisplayName', 'Average')
plot(nan, nan, 'linewidth', lw, 'color', 'black', 'linestyle', '--', 'DisplayName', 'Min/Max')


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
% COMPUTE FARM HALF WIDTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turbine = 1;
half_width_factor = 0.5;

clc; close all
figure('color', 'white')
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
    halfwidth.x(block, :) = x;
    halfwidth.massaged(block, :) = massaged_bottom_edge;
    halfwidth.raw(block, :) = raw_bottom_edge;
    
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


%% Plot Momentum thickness along x
%%% Paper figure

turbine = 1;

labelFontSize = 10;
tickFontSize = 8;

% Near, middle, and far regions
% colors(1).c = '#52b788';
% colors(2).c = '#2d6a4f';
% colors(3).c = '#081c15';

% colors(1).c = '#8390fa';
% colors(2).c = '#8390fa';
% colors(3).c = '#8390fa';

colors(1).c = '#ffbe0b';
colors(2).c = '#ffbe0b';
colors(3).c = '#ffbe0b';

data_shade_color = [0.91, 0.91 ,0.91];
% data_shade_color = '#E9E9E9';

skip = 1;
thresh_frac = 0.001;   % fraction of peak below which we consider integrand "closed"

theta = cell(1,3);
x_theta = cell(1,3);



% hold on
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
    end
end





% Plot Momentum Thickness
x_cutoff = 5.45;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5, 5, 13, 5]);
tiledlayout(1,1,'padding', 'none')
nexttile()
hold on
for block = 1:3
    x_tmp = x_theta{block};
    data_tmp = theta{block};

    data_tmp(x_tmp < x_cutoff) = nan;
    plot(x_tmp, data_tmp, 'color', colors(block).c, 'linewidth', 2, 'handlevisibility', 'off')
end
axis equal
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\theta \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% grid on

xlim([0, 52])
ylim([-6.5,2.5])

% Fit a line through blocks 2 and 3
xs = [x_theta{1}, x_theta{2}, x_theta{3}];
ys = [theta{1}, theta{2}, theta{3}];


[~, x_cutoff_idx] = min(abs(xs - x_cutoff));
xs_fit = xs(x_cutoff_idx:end);
ys_fit = ys(x_cutoff_idx:end);

P = polyfit(xs_fit, ys_fit, 1);
theta_fit = P;
x_fit = x_cutoff:50;
% x_fit = 0:50;
tmp = plot(x_fit, polyval(P, x_fit), 'color', 'black', 'linestyle', '--', 'handlevisibility', 'off');
uistack(tmp, 'bottom')

% Wake merging point
xline(x_cutoff, 'linewidth', 1, 'handlevisibility', 'off', ...
      'label', sprintf('Wake merging\n$x \\mathbin{/} D = %1.1f$', x_cutoff), 'interpreter', 'latex', ...
      'Fontsize', 6, 'LabelOrientation', 'aligned', 'LabelVerticalAlignment', 'middle', ...
      'LabelHorizontalAlignment', 'left')

% Create LaTeX string for the fit
fit_str = sprintf('$\\theta = %.3f\\,\\left(\\frac{x}{D} \\right) %+ .3f$', P(1), P(2));



% Legend
% plot(nan, nan, 'color', colors(1).c, 'linewidth', 2, 'displayname', 'Data')
% plot(nan, nan, 'color', 'k', 'linestyle', '-', 'linewidth', 1.5, 'displayname', 'Wake Merging')
plot(nan, nan, 'color', 'k', 'linestyle', '--', 'displayname', fit_str)

ax = gca;
yl = ax.YLim;

% hshade = patch([0 x_cutoff x_cutoff 0], ...
%                [yl(1) yl(1) yl(2) yl(2)], ...
%                [0.7 0.7 0.7], ...      % gray color (change if you want)
%                'FaceAlpha', 0.4, ...   % transparency
%                'EdgeColor', 'none', ...
%                'HandleVisibility', 'off');


% Shade the regions where we have data
xmins = [1, 21, 41];

for block = 1:3
    xmin = xmins(block);
    xmax = xmin + 9;

    hshade = patch([xmin xmax xmax xmin], ...
               [yl(1) yl(1) yl(2) yl(2)], ...
               data_shade_color, ...   % gray color (change if you want)
               'FaceAlpha', 1.0, ...   % transparency
               'EdgeColor', 'none', ...
               'HandleVisibility', 'off');

    % Send it behind everything
    uistack(hshade, 'bottom');
end
hold off


ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');

leg = legend('box', 'off', 'interpreter', 'latex', 'location', 'southeast', 'fontsize', 8);
leg.IconColumnWidth = 15;




% -------- DUAL TOP TICKS + BOTTOM TICKS --------
ax1 = gca;
L_over_D = 63;

% Make sure the main axis is the data axis
axes(ax1)

pause(1)
% -------------------------
% Main axis: bottom x/D + y ticks left/right + box
% -------------------------
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.Box = 'on';
ax1.TickDir = 'in';                 % bottom ticks point upward
ax1.TickLength = [0.012 0.012];
ax1.Layer = 'top';

ax1.XColor = 'k';
ax1.YColor = 'k';

% Bottom x/D ticks
ax1.XTick = 0:5:50;                 % adjust as needed

pause(1)
% -------------------------
% Axis 2: top inward x/D ticks only, no labels
% -------------------------
axTopD = axes('Position', ax1.Position, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'XLim', ax1.XLim, ...
    'YLim', ax1.YLim, ...
    'XColor', 'k', ...
    'YColor', 'none');

axTopD.XTick = ax1.XTick;
axTopD.XTickLabel = [];
axTopD.YTick = [];

axTopD.TickDir = 'in';              % top x/D ticks point downward
axTopD.TickLength = ax1.TickLength;
axTopD.LineWidth = ax1.LineWidth;
axTopD.Box = 'off';
axTopD.Layer = 'top';

pause(1)
% -------------------------
% Axis 3: top outward x/Lx ticks + labels
% -------------------------
axTopL = axes('Position', ax1.Position, ...
    'Color', 'none', ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'XLim', ax1.XLim, ...
    'YLim', ax1.YLim, ...
    'XColor', 'k', ...
    'YColor', 'none');

% Top scale values x/Lx
xtop = 0:0.1:(ax1.XLim(2)/L_over_D);
xtop_pos = xtop * L_over_D;

axTopL.XTick = xtop_pos;

labels = compose('%.1f', xtop);
labels{1} = '0';
axTopL.XTickLabel = labels;

axTopL.YTick = [];
axTopL.TickDir = 'out';             % top x/Lx ticks point upward
axTopL.TickLength = ax1.TickLength;
axTopL.LineWidth = ax1.LineWidth;
axTopL.Box = 'off';
axTopL.Layer = 'top';

% Match fonts
axTopD.FontSize = ax1.FontSize;
axTopD.FontName = ax1.FontName;
axTopD.TickLabelInterpreter = ax1.TickLabelInterpreter;

axTopL.FontSize = ax1.FontSize;
axTopL.FontName = ax1.FontName;
axTopL.TickLabelInterpreter = ax1.TickLabelInterpreter;

% Top label
hx = xlabel(axTopL, '$x \mathbin{/} L_x$', 'Interpreter', 'latex');
hx.FontSize = ax1.XLabel.FontSize;
hx.FontName = ax1.XLabel.FontName;


% Keep the original bottom label on ax1
xlabel(ax1, '$x \mathbin{/} D$', 'Interpreter', 'latex');

% Link x axes
linkaxes([ax1 axTopD axTopL], 'x');

% Return focus to main axis
axes(ax1)




% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
% fig_name = 'SingleFarm_MomentumThickness.pdf';
% pause(3)
% exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
% close all


% clear below block cutoff_idx i idx mom_integrand peak_idx search skip thresh
% clear u u1 u2 y Y turbine pk P thresh_frac tmp_data x x_cutoff c_cutoff_idx_x_fit 
% clear xs xs_fit X ys






%% Combining the (1-r)/(1+r) plots with the C_{\theta} plot


sz = 10;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,13,8]);
tiledlayout(2,1,'padding', 'loose', 'TileSpacing', 'loose')
ax1 = nexttile(1);
hold(ax1, 'on')
for block = 1:3
    x = x_theta{block};
    u1 = outer_mean(block).u;
    u2 = inner_mean(block).u;
    theta_tmp = theta{block} * 0.08;

    if block == 1
        theta_tmp(x < x_cutoff) = nan;
    end

    theta_smooth = sgolayfilt(theta_tmp, 3, 21);
    dtheta_dx = gradient(theta_smooth, x * 0.08);

    r = u2 ./ u1;
    ratio_parameter = (1 - r) ./ (1 + r);
    plot(x, ratio_parameter, 'color', colors(block).c, 'linewidth', 2, 'HandleVisibility', 'off')
end

% Legend
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', 2, 'Displayname', label)
end

hold off
box on
xlim(ax1, [0,51])
ylim(ax1, [0,0.18])
xlabel(ax1, '$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(ax1, '$\frac{1 - r}{1 + r}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Wake merging point
xline(x_cutoff, 'linewidth', 1, 'handlevisibility', 'off', ...
      'label', sprintf('Wake Merging: \n$x \\mathbin{/} D = %1.1f$', x_cutoff), 'interpreter', 'latex', ...
      'Fontsize', 8, 'LabelOrientation', 'horizontal', 'LabelVerticalAlignment', 'bottom')


yl = ax1.YLim;

% Shade the regions where we have data
xmins = [1, 21, 41];

for block = 1:3
    xmin = xmins(block);
    xmax = xmin + 9;

    hshade = patch([xmin xmax xmax xmin], ...
               [yl(1) yl(1) yl(2) yl(2)], ...
               data_shade_color, ...   % gray color (change if you want)
               'FaceAlpha', 1.0, ...   % transparency
               'EdgeColor', 'none', ...
               'HandleVisibility', 'off');

    % Send it behind everything
    uistack(hshade, 'bottom');
end
hold off

% ax = gca;
ax1.FontSize = tickFontSize;
ax1.TickLabelInterpreter = 'latex';
box(ax1, 'on')

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize);
leg.IconColumnWidth = 15;
leg.Units = 'normalized';

% 3. Move the legend using [left, bottom, width, height]
leg.Position = [0.43 0.72 0.5 0.13]; 





pause(2)

ax3 = nexttile(2);
hold(ax3, 'on')
lambda_offset = 0;
hold on
for block = 1:3
    x = x_theta{block};
    u1 = outer_mean(block).u;
    u2 = inner_mean(block).u;
    r = u2 ./ u1;
    ratio_parameter = (1 - r) ./ (1 + r);
    % theta_tmp = theta{block} * 0.08;

    % Skip if in the near-near wake
    theta_tmp = theta{block} * 0.08;
    if block == 1
        theta_tmp(x < x_cutoff) = nan;
        % lambda(x < x_cutoff) = nan;
    end

    % If there's a gap from the previous block, accumulate lambda
    % assuming r stays at the last known value across the gap
    if block > 1
        x_prev = x_theta{block-1};
        u1_prev = outer_mean(block-1).u;
        u2_prev = inner_mean(block-1).u;
        r_end = u2_prev(end) / u1_prev(end);
        r_start = r(1);
        r_gap = 0.5 * (r_end + r_start);  % average r across gap
        gap_dx = x(1) - x_prev(end);
        lambda_offset = lambda_offset + (1 - r_gap) / (1 + r_gap) * gap_dx;
    end

    % Integrate within this block, then shift by accumulated offset
    lambda = cumtrapz(x, ratio_parameter) + lambda_offset;

    % Update offset for next block
    lambda_offset = lambda(end);

    % ... rest of your plotting code
    scatter(lambda, theta_tmp, sz, 'filled', 'markerfacecolor',  colors(block).c, ...
            'markerfacealpha', 1, ...
            'HandleVisibility', 'off')


    fprintf('Min: %1.3f - Max %1.3f\n', min(theta_tmp), max(theta_tmp))
    disp(range(theta_tmp))

    % Fit a line: slope = C_theta
    P = polyfit(lambda(~isnan(theta_tmp)), theta_tmp(~isnan(theta_tmp)), 1);
    C_theta = P(1);
    
    spreading_rate(block) = C_theta;

    % Plot a line to the fit
    % x_plot = -0.05:0.01:0.6;
    perct_buff = 0.06;
    x_plot = linspace((1 - perct_buff) * min(lambda), (1 + perct_buff) * max(lambda), 5);

    if block == 1
        x_plot = linspace(0.4, 1.2, 5);
    end

    P = plot(x_plot, polyval(P, x_plot), 'linestyle', '--', 'color', colors(block).c, 'linewidth', 1.5, 'HandleVisibility', 'off');
    uistack(P, 'bottom')

    % Annotate with the fit
    if block == 1
        label = sprintf('$C_{\\theta} \\approx 0$');
        % offset = [0, -0.003];
    else
        label = sprintf('$C_{\\theta} = %.4f$', spreading_rate(block));
        % offset = [0, -0.003];
    end

    if block == 1
        x_loc = 0.6;
        y_loc = 0.024;
    elseif block == 2
        x_loc = 1.6;
        y_loc = 0.035;
    elseif block == 3
        x_loc = 2.65;
        y_loc = 0.038;
    end

    text(x_loc, y_loc, label, ...
        'Interpreter', 'latex', ...
        'FontSize', tickFontSize, ...
        'Color', colors(block).c, ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'top')

end

% % Legend
% for block = 1:3
%     label = sprintf('Block %1.0f', block);
%     plot(nan, nan, 'Color', colors(block).c, 'Linewidth', 2, 'Displayname', label)
% end
hold off



ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');

% xlabel('$\lambda$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlabel('$\lambda = \int \frac{1 - r(x)}{1 + r(x)} dx$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel('$\lambda = \displaystyle \int \frac{1 - r(x)}{1 + r(x)} \,\, dx$', ...
       'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\theta$ [m]', 'interpreter', 'latex', 'fontsize', labelFontSize)

% leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'northwest');
% leg.IconColumnWidth = 15;

box on
% axis square
xlim([0, 4])
ylim([0.015, 0.06])
yticks(0.02:0.02:0.06)
xticks(0:0.5:4)



% Add second x-axis above the top tile

% drawnow  % important: lets tiledlayout finalize axes positions
% 
% % Scaling for top axis
% scale_top = 1/63;   % example: top value = bottom value * scale_top
% 
% % Create transparent overlay axes on top of ax1
% ax2 = axes(fig, ...
%     'Position', ax1.Position, ...
%     'XAxisLocation', 'top', ...
%     'YAxisLocation', 'right', ...
%     'Color', 'none', ...
%     'YColor', 'none', ...
%     'XColor', 'k', ...
%     'Box', 'off');
% 
% % Match limits and ticks
% ax2.XLim = ax1.XLim;
% ax2.YLim = ax1.YLim;
% ax2.XTick = ax1.XTick;
% 
% % Relabel ticks with rescaled values
% xt = ax1.XTick * scale_top;
% 
% labels = strings(size(xt));
% for i = 1:numel(xt)
%     if abs(xt(i)) < 1e-12
%         labels(i) = "0";
%     else
%         labels(i) = sprintf('%.1f', xt(i));
%     end
% end
% 
% ax2.XTickLabel = labels;
% 
% % Formatting
% ax2.FontSize = tickFontSize;
% ax2.TickLabelInterpreter = 'latex';
% ax2.LineWidth = ax1.LineWidth;
% 
% xlabel(ax2, '$x \mathbin{/} L_x$', ...
%     'Interpreter', 'latex', ...
%     'FontSize', labelFontSize)
% 
% % Keep the overlay aligned if figure changes
% linkprop([ax1 ax2], {'Position', 'XLim', 'YLim'});


% Add second x-axis above the TOP tile
% Desired formatting:
%   bottom axis: x/D ticks pointing inward/upward
%   top edge:    x/D ticks pointing inward/downward
%   top axis:    x/Lx ticks + labels pointing outward/upward
%   box on with left/right ticks

drawnow  % lets tiledlayout finalize axes positions

L_over_D = 63;   % Lx / D

% Make ax1 the main data axis
axes(ax1)

% -------------------------
% Main axis: x/D on bottom, box ticks on all sides
% -------------------------
ax1.XAxisLocation = 'bottom';
ax1.YAxisLocation = 'left';
ax1.Box = 'on';
ax1.TickDir = 'in';
ax1.Layer = 'top';
ax1.LineWidth = 1;
ax1.TickLength = [0.012 0.012];

ax1.FontSize = tickFontSize;
ax1.TickLabelInterpreter = 'latex';

% Set bottom x/D ticks explicitly
ax1.XTick = 0:5:50;

% -------------------------
% Overlay axis: top x/Lx labels + outward ticks
% -------------------------
axTop = axes(fig, ...
    'Position', ax1.Position, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'Color', 'none', ...
    'YColor', 'none', ...
    'XColor', 'k', ...
    'Box', 'off', ...
    'TickDir', 'out', ...
    'TickLength', ax1.TickLength, ...
    'LineWidth', ax1.LineWidth);

% Match limits
axTop.XLim = ax1.XLim;
axTop.YLim = ax1.YLim;

% Define top-axis values in x/Lx
xtop = 0:0.1:(ax1.XLim(2) / L_over_D);

% Convert x/Lx tick values to x/D positions
xtop_pos = xtop * L_over_D;

% Keep only visible ticks
valid = xtop_pos >= ax1.XLim(1) & xtop_pos <= ax1.XLim(2);
xtop = xtop(valid);
xtop_pos = xtop_pos(valid);

% Apply top ticks
axTop.XTick = xtop_pos;

% Format labels
labels = compose('%.1f', xtop);
labels{1} = '0';
axTop.XTickLabel = labels;

% Hide y ticks from overlay axis
axTop.YTick = [];

% Match font styling
axTop.FontSize = tickFontSize;
axTop.FontName = ax1.FontName;
axTop.TickLabelInterpreter = 'latex';

% Top label
hx = xlabel(axTop, '$x \mathbin{/} L_x$', ...
    'Interpreter', 'latex', ...
    'FontSize', labelFontSize);

hx.FontName = ax1.FontName;

% Force exact alignment after tiledlayout adjustment
axTop.Units = ax1.Units;
axTop.Position = ax1.Position;
axTop.InnerPosition = ax1.InnerPosition;

% Keep overlay aligned if figure changes
hlink = linkprop([ax1 axTop], {'Position', 'InnerPosition', 'XLim', 'YLim'});
setappdata(fig, 'topAxisLink_topTile', hlink);

% Return focus to the main axis
axes(ax1)


% Save figure
save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
fig_name = 'SingleFarm_SpreadingFactor.pdf';
pause(3)
exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
close all




%% Try scaling profiles with all the computed values (all blocks overlayed)
% Optimized Tanh() + Erf()
%%% PAPER FIGURE


% Near, middle, and far regions
colors(1).c = '#fb8b24';
colors(2).c = '#d90368';
colors(3).c = '#820263';

% colors(1).c = '#6909C2';
% colors(2).c = '#B033C4';
% colors(3).c = '#C883C5';

line_alpha = 0.2;

tanh_color_trad = '#00c49a';
erf_color_trad = '#0218c2';

tanh_color_fit = '#00c49a';
erf_color_fit = '#0218c2';

legendLineWidth = 1.5;


turbine = 1;
skip = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,13,7]);
tile = tiledlayout(1,2,'padding', 'compact', 'TileSpacing', 'compact');
h(1) = nexttile();

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;
    end
end



% Legend
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', legendLineWidth, 'Displayname', label)
end

% Legend white space
plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Traditional Mixing layer fit
x_plot = -4:0.01:4;

% Tanh
plot(0.5*(1 - tanh(x_plot/2)), x_plot, 'color', tanh_color_trad, 'linewidth', 1, 'displayname', '$a = 1$', 'HandleVisibility', 'off')
plot(nan, nan, 'color', tanh_color_trad, 'linewidth', legendLineWidth, 'displayname', '$a = 1$')

% ERF
plot(0.5*(1 - erf((sqrt(pi)/2)*x_plot)), x_plot, 'color', erf_color_trad, 'linewidth', 1, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)), 'HandleVisibility', 'off')
plot(nan, nan, 'color', erf_color_trad, 'linewidth', legendLineWidth, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)))

% Inset gridlines
% Horizontal lines
P = plot([0, 1], [0, 0], 'color', 'black', 'linewidth', 1, 'HandleVisibility', 'off');
P.Color(4) = 0.5;
uistack(P, 'bottom')

P = plot([0, 1], [-2, -2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0, 1], [2, 2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Vertical lines
P = plot([0.5, 0.5], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.25, 0.25], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.75, 0.75], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Define corner coordinates (x and y)
x_square = [0 1 1 0 0];
y_square = [-4 -4 4 4 -4];

% Plot the square
bb = plot(x_square, y_square, 'black', 'LineWidth', 1.25, 'HandleVisibility', 'off');
uistack(bb, 'top')


hold off

% grid on
box on
axis square
xlim([-0.7, 1.1])
ylim([-5, 14])

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'southeast');
leg.IconColumnWidth = 15;
yticks(-8:4:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')









% Side plot zoomed in
h(2) = nexttile();
u_all   = [];
eta_all = [];

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;

        % Store only data inside the zoomed-in fit region
        fit_mask = isfinite(u_tilde) & isfinite(eta_mom) & ...
        u_tilde >= 0 & u_tilde <= 1 & ...
        eta_mom >= -4 & eta_mom <= 4;

        u_all   = [u_all; u_tilde(fit_mask)];
        eta_all = [eta_all; eta_mom(fit_mask)];
    end
end

% Fits
% Sort data by eta for cleaner plotting
[eta_all, sort_idx] = sort(eta_all);
u_all = u_all(sort_idx);

% Force fit data to double
eta_all = double(eta_all(:));
u_all   = double(u_all(:));

% Remove any remaining bad values
good = isfinite(eta_all) & isfinite(u_all);
eta_all = eta_all(good);
u_all   = u_all(good);

% -------------------------
% Linear fit: u_tilde = a eta + b
% -------------------------
P_lin = polyfit(eta_all, u_all, 1);

eta_fit = linspace(-4, 4, 300);
% u_lin_fit = polyval(P_lin, eta_fit);

a = P_lin(1);
b = P_lin(2);

% % Linear fit forced through (eta=0, u_tilde=0.5)
m = eta_all \ (u_all - 0.5);
u_lin_fit = m * eta_fit + 0.5;
% % lin_str = sprintf('$\\tilde{u} = %.3f\\,\\eta_{\\theta} + 0.5$', m);
lin_str = sprintf('$m = %.3f$', m);


lin_plot = plot(u_lin_fit, eta_fit, ...
    'color', 'black', ...
    'LineWidth', 1.5, ...
    'linestyle', '--', ...
    'DisplayName', lin_str, 'HandleVisibility', 'off');

plot(nan, nan, ...
    'color', 'black', ...
    'LineWidth', legendLineWidth, ...
    'linestyle', '--', ...
    'DisplayName', lin_str)

% -------------------------
% Tanh fit: u_tilde = 0.5*(1 - tanh(eta / (2*a)))
% -------------------------
tanh_fun = @(a, eta) 0.5*(1 - tanh(eta ./ (2*a)));

% Initial guess
a0 = 1;

% Fit
a_tanh = lsqcurvefit(tanh_fun, a0, eta_all, u_all);

u_tanh_fit = tanh_fun(a_tanh, eta_fit);

% tanh_str = sprintf('$\\frac{1}{2}\\left[1 - \\tanh\\left(\\frac{\\eta_{\\theta}}{%.3f}\\right)\\right]$', 2*a_tanh);
tanh_str = sprintf('$a = %1.3f$', a_tanh);

plot(u_tanh_fit, eta_fit, ...
    'color', tanh_color_fit, ...
    'LineWidth', 1.5, ...
    'linestyle', ':', ...
    'DisplayName', tanh_str, 'HandleVisibility', 'off')

plot(nan, nan, ...
    'color', tanh_color_fit, ...
    'LineWidth', legendLineWidth, ...
    'linestyle', ':', ...
    'DisplayName', tanh_str)


P = yline(0, 'color', 'black', 'linewidth', 1, 'alpha', 0.5, 'HandleVisibility', 'off');
uistack(P, 'bottom')


% -------------------------
% Erf fit: u_tilde = 0.5*(1 - erf(eta / sigma))
% -------------------------
erf_fun = @(sigma, eta) 0.5*(1 - erf(eta ./ sigma));

% Initial guess (canonical value is 2/sqrt(pi) ≈ 1.128)
sigma0 = 2/sqrt(pi);

% Fit
sigma_erf = lsqcurvefit(erf_fun, sigma0, eta_all, u_all);

u_erf_fit = erf_fun(sigma_erf, eta_fit);

% erf_str = sprintf('$\\frac{1}{2}\\left[1 - \\mathrm{erf}\\left(\\frac{\\eta_{\\theta}}{%.3f}\\right)\\right]$', sigma_erf);
erf_str = sprintf('$\\sigma = %1.3f$', sigma_erf);

plot(u_erf_fit, eta_fit, ...
    'color', erf_color_fit, ...
    'LineWidth', 1.5, ...
    'linestyle', ':', ...
    'DisplayName', erf_str, 'HandleVisibility', 'off')

plot(nan, nan, ...
    'color', erf_color_fit, ...
    'LineWidth', legendLineWidth, ...
    'linestyle', ':', ...
    'DisplayName', erf_str)


uistack(lin_plot, 'top')
hold off
axis square
xlim([0, 1])
ylim([-4,4])
xticks(0:0.25:1)

uistack(findobj(ax, 'Type', 'rectangle'), 'top'); % keep border visible

box on
grid on
leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'northwest');
leg.IconColumnWidth = 15;
yticks(-8:2:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')


xlabel(tile, '$\tilde{u} = \frac{u - u_{2}}{u_{1} - u_{2}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$\eta_{\theta} = \frac{z - \delta_{1/2}}{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
% fig_name = 'SingleFarm_MixingLayerScaling.pdf';
% pause(3)
% exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
% close all






%% Try scaling profiles with all the computed values (all blocks overlayed)
% Optimized Tanh() + Erf() Optimized per block
%%% PAPER FIGURE


% Near, middle, and far regions
colors(1).c = '#fb8b24';
colors(2).c = '#d90368';
colors(3).c = '#820263';

% colors(1).c = '#6909C2';
% colors(2).c = '#B033C4';
% colors(3).c = '#C883C5';

line_alpha = 0.2;

tanh_color_trad = '#00c49a';
erf_color_trad = '#0218c2';

tanh_color_fit = '#00c49a';
erf_color_fit = '#0218c2';

legendLineWidth = 1.5;


turbine = 1;
skip = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,13,7]);
tile = tiledlayout(1,2,'padding', 'compact', 'TileSpacing', 'compact');
h(1) = nexttile();

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;
    end
end



% Legend
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', legendLineWidth, 'Displayname', label)
end

% Legend white space
plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Traditional Mixing layer fit
x_plot = -4:0.01:4;

% Tanh
plot(0.5*(1 - tanh(x_plot/2)), x_plot, 'color', tanh_color_trad, 'linewidth', 1, 'displayname', '$a = 1$', 'HandleVisibility', 'off')
plot(nan, nan, 'color', tanh_color_trad, 'linewidth', legendLineWidth, 'displayname', '$a = 1$')

% ERF
plot(0.5*(1 - erf((sqrt(pi)/2)*x_plot)), x_plot, 'color', erf_color_trad, 'linewidth', 1, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)), 'HandleVisibility', 'off')
plot(nan, nan, 'color', erf_color_trad, 'linewidth', legendLineWidth, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)))

% Inset gridlines
% Horizontal lines
P = plot([0, 1], [0, 0], 'color', 'black', 'linewidth', 1, 'HandleVisibility', 'off');
P.Color(4) = 0.5;
uistack(P, 'bottom')

P = plot([0, 1], [-2, -2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0, 1], [2, 2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Vertical lines
P = plot([0.5, 0.5], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.25, 0.25], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.75, 0.75], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Define corner coordinates (x and y)
x_square = [0 1 1 0 0];
y_square = [-4 -4 4 4 -4];

% Plot the square
bb = plot(x_square, y_square, 'black', 'LineWidth', 1.25, 'HandleVisibility', 'off');
uistack(bb, 'top')


hold off

% grid on
box on
axis square
xlim([-0.7, 1.1])
ylim([-5, 14])

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'southeast');
leg.IconColumnWidth = 15;
yticks(-8:4:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')









% Side plot zoomed in
% Side plot zoomed in
h(2) = nexttile();

% Storage for per-block fit parameters
a_tanh_block = nan(1,3);
sigma_erf_block = nan(1,3);
m_block = nan(1,3);

hold on
for block = 1:3

    u_block   = [];
    eta_block = [];

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;

        % Store only data inside the zoomed-in fit region
        fit_mask = isfinite(u_tilde) & isfinite(eta_mom) & ...
            u_tilde >= 0 & u_tilde <= 1 & ...
            eta_mom >= -4 & eta_mom <= 4;

        u_block   = [u_block; u_tilde(fit_mask)];
        eta_block = [eta_block; eta_mom(fit_mask)];
    end

    % Clean and sort
    eta_block = double(eta_block(:));
    u_block   = double(u_block(:));
    good = isfinite(eta_block) & isfinite(u_block);
    eta_block = eta_block(good);
    u_block   = u_block(good);
    [eta_block, sort_idx] = sort(eta_block);
    u_block = u_block(sort_idx);

    eta_fit = linspace(-4, 4, 300);

    % ---- Linear fit forced through (eta=0, u_tilde=0.5) ----
    m_block(block) = eta_block \ (u_block - 0.5);

    % ---- Tanh fit ----
    tanh_fun = @(a, eta) 0.5*(1 - tanh(eta ./ (2*a)));
    a_tanh_block(block) = lsqcurvefit(tanh_fun, 1, eta_block, u_block, [], [], ...
        optimoptions('lsqcurvefit','Display','off'));

    % ---- Erf fit ----
    erf_fun = @(sigma, eta) 0.5*(1 - erf(eta ./ sigma));
    sigma_erf_block(block) = lsqcurvefit(erf_fun, 2/sqrt(pi), eta_block, u_block, [], [], ...
        optimoptions('lsqcurvefit','Display','off'));

    % Plot fits with block color, different line styles
    u_tanh_fit = tanh_fun(a_tanh_block(block), eta_fit);
    u_erf_fit  = erf_fun(sigma_erf_block(block), eta_fit);
    u_lin_fit  = m_block(block) * eta_fit + 0.5;

    plot(u_tanh_fit, eta_fit, 'color', colors(block).c, 'LineWidth', 1.5, ...
        'linestyle', '-', 'HandleVisibility', 'off')
    plot(u_erf_fit, eta_fit, 'color', colors(block).c, 'LineWidth', 1.5, ...
        'linestyle', ':', 'HandleVisibility', 'off')
    plot(u_lin_fit, eta_fit, 'color', colors(block).c, 'LineWidth', 1.5, ...
        'linestyle', '--', 'HandleVisibility', 'off')
end

% Centerline
P = yline(0, 'color', 'black', 'linewidth', 1, 'alpha', 0.5, 'HandleVisibility', 'off');
uistack(P, 'bottom')

% ---- Legend entries ----
% Block colors
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', legendLineWidth, 'Displayname', label)
end
plot(nan, nan, 'color', 'white', 'displayname', ' ')  % spacer

% Line style legend
plot(nan, nan, 'color', 'k', 'LineWidth', legendLineWidth, 'linestyle', '-',  'DisplayName', 'tanh fit')
plot(nan, nan, 'color', 'k', 'LineWidth', legendLineWidth, 'linestyle', ':',  'DisplayName', 'erf fit')
plot(nan, nan, 'color', 'k', 'LineWidth', legendLineWidth, 'linestyle', '--', 'DisplayName', 'linear fit')

hold off
axis square
xlim([0, 1])
ylim([-4, 4])
xticks(0:0.25:1)
box on
grid on
leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'northwest');
leg.IconColumnWidth = 15;
yticks(-8:2:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')

% Print fit parameters to command window
fprintf('\n--- Per-block fit parameters ---\n')
for block = 1:3
    fprintf('Block %d:  a_tanh = %.4f,  sigma_erf = %.4f,  m_lin = %.4f\n', ...
        block, a_tanh_block(block), sigma_erf_block(block), m_block(block))
end



% uistack(lin_plot, 'top')
hold off
axis square
xlim([0, 1])
ylim([-4,4])
xticks(0:0.25:1)

uistack(findobj(ax, 'Type', 'rectangle'), 'top'); % keep border visible

box on
grid on
leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'northwest');
leg.IconColumnWidth = 15;
yticks(-8:2:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')


xlabel(tile, '$\tilde{u} = \frac{u - u_{2}}{u_{1} - u_{2}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$\eta_{\theta} = \frac{z - \delta_{1/2}}{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
% fig_name = 'SingleFarm_MixingLayerScaling.pdf';
% pause(3)
% exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
% close all




%% Try scaling profiles with all the computed values (all blocks overlayed)
% Piecewise linear fits
%%% PAPER FIGURE


% Near, middle, and far regions
colors(1).c = '#fb8b24';
colors(2).c = '#d90368';
colors(3).c = '#820263';

% colors(1).c = '#6909C2';
% colors(2).c = '#B033C4';
% colors(3).c = '#C883C5';

line_alpha = 0.25;

tanh_color_trad = '#00c49a';
erf_color_trad = '#0218c2';

% tanh_color_fit = '#ffbe0a';
% erf_color_fit = '#ce95db';
tanh_color_fit = '#00c49a';
erf_color_fit = '#0218c2';

legendLineWidth = 1.5;


turbine = 1;
skip = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,13,7]);
tile = tiledlayout(1,2,'padding', 'compact', 'TileSpacing', 'compact');
h(1) = nexttile();

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;
    end
end



% Legend
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', legendLineWidth, 'Displayname', label)
end

% Legend white space
plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Traditional Mixing layer fit
x_plot = -4:0.01:4;
% Tanh
% plot(0.5*(1 - tanh(x_plot/2)), x_plot, 'color', 'blue', 'linewidth', 1, 'displayname', '$\frac{1}{2} \left[ 1 - \tanh \left( \frac{\eta_{\theta}}{2} \right) \right]$')
plot(0.5*(1 - tanh(x_plot/2)), x_plot, 'color', tanh_color_trad, 'linewidth', 1, 'displayname', '$a = 1$', 'HandleVisibility', 'off')
plot(nan, nan, 'color', tanh_color_trad, 'linewidth', legendLineWidth, 'displayname', '$a = 1$')

% ERF
% plot(0.5*(1 - erf((sqrt(pi)/2)*x_plot)), x_plot, 'color', 'red', 'linewidth', 1, 'DisplayName', '$\frac{1}{2} \left[ 1 - erf \left( \frac{\sqrt{\pi}\eta_{\theta}}{2} \right) \right]$')
plot(0.5*(1 - erf((sqrt(pi)/2)*x_plot)), x_plot, 'color', erf_color_trad, 'linewidth', 1, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)), 'HandleVisibility', 'off')
plot(nan, nan, 'color', erf_color_trad, 'linewidth', legendLineWidth, 'DisplayName', sprintf('$\\sigma = %1.3f$', 2 / sqrt(pi)))

% Inset gridlines
% Horizontal lines
P = plot([0, 1], [0, 0], 'color', 'black', 'linewidth', 1, 'HandleVisibility', 'off');
P.Color(4) = 0.5;
uistack(P, 'bottom')

P = plot([0, 1], [-2, -2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0, 1], [2, 2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Vertical lines
P = plot([0.5, 0.5], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.25, 0.25], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

P = plot([0.75, 0.75], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
P.Color(4) = 0.25;
uistack(P, 'bottom')

% Define corner coordinates (x and y)
x_square = [0 1 1 0 0];
y_square = [-4 -4 4 4 -4];

% Plot the square
bb = plot(x_square, y_square, 'black', 'LineWidth', 1.25, 'HandleVisibility', 'off');
uistack(bb, 'top')


hold off

% grid on
box on
axis square
xlim([-0.7, 1.1])
ylim([-5, 14])

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'southeast');
leg.IconColumnWidth = 15;
yticks(-8:4:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')









% Side plot zoomed in
h(2) = nexttile();
u_all   = [];
eta_all = [];

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    u(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

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
        P = plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;

        % Store only data inside the zoomed-in fit region
        fit_mask = isfinite(u_tilde) & isfinite(eta_mom) & ...
        u_tilde >= 0 & u_tilde <= 1 & ...
        eta_mom >= -4 & eta_mom <= 4;

        u_all   = [u_all; u_tilde(fit_mask)];
        eta_all = [eta_all; eta_mom(fit_mask)];
    end
end

% Plot the square
% plot(x_square, y_square, 'black', 'LineWidth', 1, 'HandleVisibility', 'off')

% Fits
% Sort data by eta for cleaner plotting
[eta_all, sort_idx] = sort(eta_all);
u_all = u_all(sort_idx);

% Force fit data to double
eta_all = double(eta_all(:));
u_all   = double(u_all(:));

% Remove any remaining bad values
good = isfinite(eta_all) & isfinite(u_all);
eta_all = eta_all(good);
u_all   = u_all(good);

% -------------------------
% Linear fit: u_tilde = a eta + b
% -------------------------
P_lin = polyfit(eta_all, u_all, 1);

eta_fit = linspace(-4, 4, 300);
% u_lin_fit = polyval(P_lin, eta_fit);

a = P_lin(1);
b = P_lin(2);

% % Linear fit forced through (eta=0, u_tilde=0.5)
m = eta_all \ (u_all - 0.5);
u_lin_fit = m * eta_fit + 0.5;
% % lin_str = sprintf('$\\tilde{u} = %.3f\\,\\eta_{\\theta} + 0.5$', m);
lin_str = sprintf('$m = %.3f$', m);


lin_plot = plot(u_lin_fit, eta_fit, ...
    'color', 'black', ...
    'LineWidth', 1.5, ...
    'linestyle', '--', ...
    'DisplayName', lin_str, 'HandleVisibility', 'off');

plot(nan, nan, ...
    'color', 'black', ...
    'LineWidth', legendLineWidth, ...
    'linestyle', '--', ...
    'DisplayName', lin_str)



%%% Piece-wise linear
% Split data into three regions
% outer = u_all >= 0.75;
% core  = u_all >= 0.25 & u_all < 0.55;
% inner = u_all < 0.25;

% Classify by eta
outer = eta_all <= -1;
core = eta_all >= -1 & eta_all <= 1;
inner = eta_all >= 1;


m_outer = polyfit(eta_all(outer), u_all(outer), 1);
m_core  = polyfit(eta_all(core),  u_all(core),  1);
m_inner = polyfit(eta_all(inner), u_all(inner), 1);


% Colors based on region
inner_color = '#fee440';
core_color = '#00bbf9';
outer_color = '#00f5d4';

plot(nan, nan, 'color', 'white', 'displayname', ' ')
% u = 0 - 0.25
eta_tmp = linspace(1.25, 4, 10);
plot(eta_tmp * m_inner(1) + m_inner(2), eta_tmp, 'color', inner_color, 'linewidth', 2, ...
     'displayname', sprintf('$m = %1.3f$', m_inner(1)))

% u = 0.25 - 0.75
eta_tmp = linspace(-0.75, 0.75, 10);
plot(eta_tmp * m_core(1) + m_core(2), eta_tmp, 'color', core_color, 'linewidth', 2, ...
     'displayname', sprintf('$m = %1.3f$', m_core(1)))

% u = 0.75 - 1
eta_tmp = linspace(-4, -1.25, 10);
plot(eta_tmp * m_outer(1) + m_outer(2), eta_tmp, 'color', outer_color, 'linewidth', 2, ...
     'displayname', sprintf('$m = %1.3f$', m_outer(1)))


% Print slopes from different regions
fprintf('\n\nInner slope: %1.3f, Core slope: %1.3f, Outer slope: %1.3f\n', m_inner(1), m_core(1), m_outer(1))
fprintf('Total fit slope: %1.3f\n\n', m)

% yline([-1, 1], 'handlevisibility', 'off')


uistack(lin_plot, 'top')
hold off
axis square
xlim([0, 1])
ylim([-4,4])
xticks(0:0.25:1)

uistack(findobj(ax, 'Type', 'rectangle'), 'top'); % keep border visible

box on
grid on
leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'northwest');
leg.IconColumnWidth = 15;
yticks(-8:2:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')


xlabel(tile, '$\tilde{u} = \frac{u - u_{2}}{u_{1} - u_{2}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$\eta_{\theta} = \frac{z - \delta_{1/2}}{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
% fig_name = 'SingleFarm_MixingLayerScaling.pdf';
% pause(3)
% exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
% close all




%% Try scaling profiles with all the computed values (all blocks overlayed)
% MIXING LAYER SCALED SHEAR STRESS (u'w')
%%% PAPER FIGURE


% Near, middle, and far regions
colors(1).c = '#fb8b24';
colors(2).c = '#d90368';
colors(3).c = '#820263';

% colors(1).c = '#6909C2';
% colors(2).c = '#B033C4';
% colors(3).c = '#C883C5';

line_alpha = 0.5;

tanh_color_trad = '#00c49a';
erf_color_trad = '#0218c2';

tanh_color_fit = '#00c49a';
erf_color_fit = '#0218c2';

legendLineWidth = 1.5;


turbine = 1;
skip = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [5,5,13,7]);
tile = tiledlayout(1,2,'padding', 'compact', 'TileSpacing', 'compact');
h(1) = nexttile();

hold on
for block = 1:3

    % Velocity data
    u = cleaned(block).u;
    uv = data(block).uv;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    uv(Y > 3) = nan;

    % Scale different profiles
    for i = 1:skip:195

        % Skip if in the near-near wake
        if block == 1 && x(i) < x_cutoff
            continue
        end

        uv_profile = uv(:, i);
        half_width_local = half_width(i);
        momentum_thickness_local = momentum_thickness(i);

        % Local inner/outer flows
        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);

        % Scaled coordinate
        eta_mom = (y - half_width_local) / momentum_thickness_local;

        % Scaled velocity
        % uv_tilde = uv_profile ./ ((u1 - u2).^2);
        uv_tilde = uv_profile/ (u1^2);

        P = plot(uv_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off');
        P.Color(4) = line_alpha;
    end
end



% Legend
for block = 1:3
    label = sprintf('Block %1.0f', block);
    plot(nan, nan, 'Color', colors(block).c, 'Linewidth', legendLineWidth, 'Displayname', label)
end

% % Inset gridlines
% % Horizontal lines
% P = plot([0, 1], [0, 0], 'color', 'black', 'linewidth', 1, 'HandleVisibility', 'off');
% P.Color(4) = 0.5;
% uistack(P, 'bottom')
% 
% P = plot([0, 1], [-2, -2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
% P.Color(4) = 0.25;
% uistack(P, 'bottom')
% 
% P = plot([0, 1], [2, 2], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
% P.Color(4) = 0.25;
% uistack(P, 'bottom')
% 
% % Vertical lines
% P = plot([0.5, 0.5], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
% P.Color(4) = 0.25;
% uistack(P, 'bottom')
% 
% P = plot([0.25, 0.25], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
% P.Color(4) = 0.25;
% uistack(P, 'bottom')
% 
% P = plot([0.75, 0.75], [-4, 4], 'color', 'black', 'linewidth', 0.25, 'HandleVisibility', 'off');
% P.Color(4) = 0.25;
% uistack(P, 'bottom')

% Define corner coordinates (x and y)
% x_square = [0 1 1 0 0];
% y_square = [-4 -4 4 4 -4];

% x_square = [-0.25 0.5 0.5 -0.25 -0.25];
% y_square = [-4 -4 6 6 -4];
% 
% % Plot the square
% bb = plot(x_square, y_square, 'black', 'LineWidth', 1.25, 'HandleVisibility', 'off');
% uistack(bb, 'top')



box on
axis square
xlim([-10E-3, 10E-3])
ylim([-4, 14])

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize, 'location', 'southwest');
leg.IconColumnWidth = 15;
yticks(-8:4:16)

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')

xlabel("$\overline{u'w'} \mathbin{/} u_{\infty}^2$", 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\eta_{\theta} = \frac{z - \delta_{1/2}}{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)






% Side plot of eddy viscosity
h(2) = nexttile();

% Common eta grid for averaging
eta_common = linspace(-4, 6, 200)';
uv_collected = [];

hold on
for block = 2:3

    % Velocity data
    u = cleaned(block).u;
    uv = data(block).uv;
    half_width = halfwidth.massaged(block,:);
    momentum_thickness = theta{block};
    x = x_theta{block};

    X = data(block).X;
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Bulk crop data to focus on mixing layer
    uv(Y > 3) = nan;

    % Scale different profiles
    for i = 1:size(uv, 2)

        uv_profile = uv(:, i);
        half_width_local = half_width(i);
        momentum_thickness_local = momentum_thickness(i);

        % Local inner/outer flows
        u1 = outer_mean(block).u(i);
        u2 = inner_mean(block).u(i);

        % Scaled coordinate
        eta_mom = (y - half_width_local) / momentum_thickness_local;

        % Scaled stress (not normalized)
        uv_tilde = uv_profile / (u1^2);

        % Interpolate onto common eta grid
        good = isfinite(uv_tilde) & isfinite(eta_mom);
        if nnz(good) >= 3
            uv_interp = interp1(eta_mom(good), uv_tilde(good), eta_common, 'linear', NaN);
            uv_collected = [uv_collected, uv_interp];
        end
    end
end

% Average across all profiles
uv_mean = mean(uv_collected, 2, 'omitnan');
uv_std  = std(uv_collected, 0, 2, 'omitnan');

% Plot mean with shaded uncertainty
fill(h(1), [uv_mean - uv_std; flipud(uv_mean + uv_std)], ...
     [eta_common; flipud(eta_common)], ...
     [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.4, ...
     'HandleVisibility', 'off')

plot(h(1), uv_mean, eta_common, 'color', 'black', 'linewidth', 1, ...
    'DisplayName', 'Mean')
% $\langle u''w'' \rangle / u_{\infty}^2$
hold off



% Compute eddy viscosity
% Collect mean quantities
U1_mean    = mean([outer_mean(2).u, outer_mean(3).u], 'omitnan');
U2_mean    = mean([inner_mean(2).u, inner_mean(3).u], 'omitnan');
DU_mean    = U1_mean - U2_mean;
theta_mean = mean([theta{2}, theta{3}], 'omitnan') * 0.08;  % meters
stress_mean = mean(uv_mean(eta_common >= -4 & eta_common <= 4), 'omitnan');

fprintf('=== Inputs ===\n')
fprintf('U1      = %.3f m/s\n', U1_mean)
fprintf('U2      = %.3f m/s\n', U2_mean)
fprintf('DU      = %.3f m/s\n', DU_mean)
fprintf('theta   = %.5f m\n', theta_mean)
fprintf('stress  = %.6f (normalized by U1^2)\n', stress_mean)
fprintf('m_outer = %.4f\n', m_outer(1))
fprintf('m_core  = %.4f\n', m_core(1))
fprintf('m_inner = %.4f\n', m_inner(1))

fprintf('\n=== Eddy Viscosity ===\n')
nu_T_outer = -stress_mean * U1_mean^2 * theta_mean / (DU_mean * m_outer(1));
nu_T_core  = -stress_mean * U1_mean^2 * theta_mean / (DU_mean * m_core(1));
nu_T_inner = -stress_mean * U1_mean^2 * theta_mean / (DU_mean * m_inner(1));

fprintf('nu_T outer: %.4e m^2/s\n', nu_T_outer)
fprintf('nu_T core:  %.4e m^2/s\n', nu_T_core)
fprintf('nu_T inner: %.4e m^2/s\n', nu_T_inner)

fprintf('\n=== Ratios ===\n')
fprintf('outer/core: %.2f\n', nu_T_outer / nu_T_core)
fprintf('inner/core: %.2f\n', nu_T_inner / nu_T_core)
fprintf('inner/outer: %.2f\n', nu_T_inner / nu_T_outer)

fprintf('\n=== Non-dimensional (nu_T / U1*theta) ===\n')
fprintf('outer: %.4f\n', nu_T_outer / (U1_mean * theta_mean))
fprintf('core:  %.4f\n', nu_T_core  / (U1_mean * theta_mean))
fprintf('inner: %.4f\n', nu_T_inner / (U1_mean * theta_mean))



% Plot eddy viscosities as a barchart
hold on


nu_T_vals = [nu_T_outer, nu_T_core, nu_T_inner];
labels = {'Outer', 'Mixing', 'Inner'};
% bar_colors = [0 0 1; 0 0 0; 1 0 0];  % blue, black, red — match piecewise fit lines
bar_colors = [hex2rgb(outer_color); hex2rgb(core_color); hex2rgb(inner_color)];

b = barh(1:3, nu_T_vals, 0.5, 'FaceColor', 'flat', 'CData', bar_colors, 'EdgeColor', 'none');
b.FaceColor = 'flat';
% for i = 1:3
%     b.CData(i,:) = hex2rgb(colors(i).c);
% end

set(gca, 'YTick', 1:3, 'YTickLabel', labels)
set(gca, 'YDir', 'reverse')  % outer on top, inner on bottom — matches spatial layout
xlabel('$\nu_T$ [m$^2$/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Add value labels at end of each bar
for k = 1:3
    text(nu_T_vals(k) + 0.001, k, sprintf('%.3f', nu_T_vals(k)), ...
        'VerticalAlignment', 'middle', ...
        'FontSize', tickFontSize, ...
        'Interpreter', 'latex')
end

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex')
box on

ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex')

hold off
axis square
xlim([0, 0.05])
ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');
set(gca, 'YDir', 'reverse')



% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MixingLayer';
% fig_name = 'SingleFarm_MixingLayerScaling.pdf';
% pause(3)
% exportgraphics(fig, fullfile(save_folder, fig_name), 'resolution', 600)
% close all




%% Try scaling profiles with all the computed values (all blocks overlayed) CALUDE BOI
% Save all valid profiles and fit in -4 <= eta_theta <= 4

% turbine = 1;
% skip = 5;
% 
% clc; close all
% figure('color', 'white', 'units', 'centimeters', 'position', [5,5,10,8])
% tiledlayout(1,1,'padding', 'none')
% nexttile()
% hold on
% 
% % Containers for all profiles
% all_eta = [];
% all_utilde = [];
% 
% for block = 1:3
% 
%     % Velocity data
%     u = cleaned(block).u;
%     half_width = halfwidth.massaged(block,:);
%     momentum_thickness = theta{block};
%     x = x_theta{block};
% 
%     X = data(block).X; %#ok<NASGU>
%     Y = data(block).Y - turbine_positions(turbine);
%     y = Y(:,1);
% 
%     % Bulk crop data to focus on mixing layer
%     u(Y > 3) = nan;
% 
%     % Safer number of slices
%     nSlices = size(u,2);
% 
%     % Scale different profiles
%     for i = 1:skip:nSlices
% 
%         % Skip if in the near-near wake
%         if block == 1 && x(i) < x_cutoff
%             continue
%         end
% 
%         u_profile = u(:, i);
%         half_width_local = half_width(i);
%         momentum_thickness_local = momentum_thickness(i);
% 
%         % Skip bad scaling values
%         if isnan(half_width_local) || isnan(momentum_thickness_local) || momentum_thickness_local == 0
%             continue
%         end
% 
%         % Local inner/outer flows
%         u1 = outer_mean(block).u(i);
%         u2 = inner_mean(block).u(i);
% 
%         if isnan(u1) || isnan(u2) || (u1 - u2) == 0
%             continue
%         end
% 
%         % Scaled coordinate
%         eta_mom = (y - half_width_local) / momentum_thickness_local;
% 
%         % Scaled velocity
%         u_tilde = (u_profile - u2) / (u1 - u2);
% 
%         % Keep only finite values for this profile
%         good = isfinite(eta_mom) & isfinite(u_tilde);
%         eta_mom = eta_mom(good);
%         u_tilde = u_tilde(good);
% 
%         % Plot profile
%         plot(u_tilde, eta_mom, 'color', colors(block).c, 'HandleVisibility', 'off')
% 
%         % Save all points
%         all_eta = [all_eta; eta_mom(:)];
%         all_utilde = [all_utilde; u_tilde(:)];
%     end
% end
% 
% % Restrict to fit window: -4 <= eta_theta <= 4
% fit_mask = isfinite(all_eta) & isfinite(all_utilde) & (all_eta >= -4) & (all_eta <= 4);
% 
% eta_fit = all_eta(fit_mask);
% u_fit = all_utilde(fit_mask);
% 
% % Optional: sort by eta for plotting smooth curves
% [eta_fit_sorted, sort_idx] = sort(eta_fit);
% u_fit_sorted = u_fit(sort_idx);
% 
% % Define fit models: u_tilde = f(eta)
% tanh_model = @(p,eta) 0.5 * (1 - tanh(p(1)*eta + p(2)));
% erf_model  = @(p,eta) 0.5 * (1 - erf( p(1)*eta + p(2)));
% 
% % Objective functions
% obj_tanh = @(p) sum((u_fit - tanh_model(p, eta_fit)).^2);
% obj_erf  = @(p) sum((u_fit - erf_model(p, eta_fit)).^2);
% 
% % Initial guesses
% p0_tanh = [1, 0];
% p0_erf  = [sqrt(pi)/2, 0];
% 
% % Nonlinear fits using fminsearch
% opts = optimset('Display','off');
% p_tanh = fminsearch(obj_tanh, p0_tanh, opts);
% p_erf  = fminsearch(obj_erf,  p0_erf,  opts);
% 
% % Linear fit: u_tilde = m*eta + b
% P_lin = polyfit(eta_fit, u_fit, 1);
% 
% % Evaluate fits on a smooth eta grid
% eta_plot = linspace(-4, 4, 400);
% 
% u_tanh = tanh_model(p_tanh, eta_plot);
% u_erf  = erf_model(p_erf, eta_plot);
% u_lin  = polyval(P_lin, eta_plot);
% 
% % Plot the three fits
% % h1 = plot(u_tanh, eta_plot, 'magenta',  'LineWidth', 1.5, 'DisplayName', 'tanh fit');
% % h2 = plot(u_erf,  eta_plot, 'red', 'LineWidth', 1.5, 'DisplayName', 'erf fit');
% % h3 = plot(u_lin,  eta_plot, 'cyan',  'LineWidth', 1.8, 'DisplayName', 'linear fit');
% 
% % uistack([h1 h2 h3], 'top')
% 
% % Define corner coordinates (x and y)
% x_square = [0 1 1 0 0];
% y_square = [-4 -4 4 4 -4];
% 
% % Plot the square
% plot(x_square, y_square, 'black', 'LineWidth', 1, 'HandleVisibility', 'off')
% 
% % Legend entries for blocks
% for block = 1:3
%     label = sprintf('Block %1.0f', block);
%     plot(nan, nan, 'Color', colors(block).c, 'Linewidth', 2, 'Displayname', label)
% end
% 
% % Legend white space
% plot(nan, nan, 'color', 'white', 'displayname', ' ')
% 
% axis square
% xlim([-1, 1.1])
% box on
% grid on
% leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize);
% leg.Layout.Tile = 'east';
% yticks(-8:4:16)
% 
% ax = gca;
% ax.FontSize = tickFontSize;
% set(gca, 'TickLabelInterpreter', 'latex');
% set(gca, 'YDir', 'reverse')
% 
% xlabel('$\tilde{u}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$\eta_{\theta}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% 
% % Display fitted equations in command window
% fprintf('\nTanh fit:\n')
% fprintf('u_tilde = 0.5 * (1 - tanh(%.4f * eta_theta %+ .4f))\n', p_tanh(1), p_tanh(2))
% 
% fprintf('\nErf fit:\n')
% fprintf('u_tilde = 0.5 * (1 - erf(%.4f * eta_theta %+ .4f))\n', p_erf(1), p_erf(2))
% 
% fprintf('\nLinear fit:\n')
% fprintf('u_tilde = %.4f * eta_theta %+ .4f\n', P_lin(1), P_lin(2))







