% Attempt to show the far wake can be scaled: using mixing layer approach
% Zein Sadek

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
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
% experiment = {'Farm2Farm_10D_Gap'};
blocks = [1,2,3];

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
% CREATE A MASSAGED VERSION OF THE DATA AT THE SEAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

for block = 1:3
    for c = 1:length(components)
        component = components{c};

        % Pull in block data
        u = data(block).(component);
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

clear components block c component u X Y blended left right discontinuity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SEGMENT TURBINES, GET CENTERLINE VELOCITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turbine/Farm dimensions
D = 80;
Sy = 3;
u_inf = 7.5;

% T1 ~ edge turbine
% T4 ~ center turbine
turbine_positions = -9:3:3;

% Loop through and crop different wakes
for t = 1:length(turbine_positions)
    for block = 1:3

        % Load data
        X = data(block).X;
        Y = data(block).Y;
        x = X(1,:);
        y = Y(:,1);
        
        raw_u = data(block).u;
        massaged_u = cleaned(block).u;

        centered_y = y - turbine_positions(t);
        mask = (centered_y >= -Sy) & (centered_y <= Sy);
    
        % Save
        wakes(t, block).raw.u = raw_u(mask, :);
        wakes(t, block).massaged.u = massaged_u(mask, :);
        wakes(t, block).Y = Y(mask, :);
        wakes(t, block).X = X(mask, :);

        % Raw + Massaged
        [r,~] = size(wakes(t, block).raw.u);
        center_index = round(r/2);
    
        % Centerline of both datasets
        raw_center_velocity =  wakes(t, block).raw.u(center_index, :);
        massaged_center_velocity =  wakes(t, block).massaged.u(center_index, :);

        wakes(t, block).raw.centerline_velocity = raw_center_velocity;
        wakes(t, block).massaged.centerline_velocity = massaged_center_velocity;
    end
end

% clear t centered_y mask block X Y y raw_u massaged_u centered_y mask
% clear center_index massaged_center_velocity r raw_center_velocity




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CENTERLINE VELOCITIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear mean_centerline_velocity

% Turbine colors
turbine_colors = slanCM('imola', 5);

% Plot centerlines 
clc; close all
figure('color', 'white')
tiledlayout(2,1)

h(1) = nexttile;
hold on
for block = 1:3
    % Average over all turbines per block
    tmp = nan(5, length(wakes(1,1).X));
    for t = 2:5
        plot(wakes(t, block).X(1,:), wakes(t, block).raw.centerline_velocity, 'color', turbine_colors(t,:))
        tmp(t,:) = wakes(t, block).raw.centerline_velocity;
    end
    plot(wakes(t, block).X(1,:), mean(tmp, 1, 'omitnan'), 'color', 'black', 'linewidth', 2)

    % Save mean
    mean_centerline_velocity(block).raw.u = mean(tmp, 1, 'omitnan');
    mean_centerline_velocity(block).raw.x = wakes(t, block).X(1,:);
end
hold off
title('Raw')


h(1) = nexttile;
hold on
for block = 1:3
    % Average over all turbines per block
    tmp = nan(5, length(wakes(1,1).X));
    for t = 2:5
        plot(wakes(t, block).X(1,:), wakes(t, block).massaged.centerline_velocity, 'color', turbine_colors(t,:))
        tmp(t,:) = wakes(t, block).massaged.centerline_velocity;
    end
    plot(wakes(t, block).X(1,:), mean(tmp, 1, 'omitnan'), 'color', 'black', 'linewidth', 2)

    % Save mean
    mean_centerline_velocity(block).massaged.u = mean(tmp, 1, 'omitnan');
    mean_centerline_velocity(block).massaged.x = wakes(t, block).X(1,:);
end
hold off
title('Massaged')


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE FARM HALF WIDTH
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

u_inf = 7.5;
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

    % Get centerline velocities
    raw_center_velocity = mean_centerline_velocity(block).raw.u;
    massaged_center_velocity = mean_centerline_velocity(block).massaged.u;

    % Find half-width velocity
    raw_center_half_velocity = half_width_factor * (raw_center_velocity + u_inf);
    massaged_center_half_velocity = half_width_factor * (massaged_center_velocity + u_inf);

    % Inside half-width
    raw_wake(raw_wake > raw_center_half_velocity) = nan;
    massaged_wake(massaged_wake > massaged_center_half_velocity) = nan;

    % Detect top and bottom half-widths
    raw_isnan = isnan(raw_wake);
    massaged_isnan = isnan(massaged_wake);



    %%% Top (bottom unflipped)
    % Y_bottom = wakes(t).Y(center_index:end, 1);
    raw_isnan_bottom = raw_isnan;
    massaged_isnan_bottom = massaged_isnan;

    % Vertical difference to find edge
    raw_isnan_bottom_diff = abs(diff(raw_isnan_bottom, 1, 1));
    massaged_isnan_bottom_diff = abs(diff(massaged_isnan_bottom, 1, 1));
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);
    raw_bottom_edge = Y(raw_bottom_edge_index);
    massaged_bottom_edge = Y(massaged_bottom_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y(1);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y(1);

    % u_tilde = (massaged_wake_uncropped - massaged_center_velocity) ./ (10 - massaged_center_velocity);

    % Plot
    contourf(X, Y, massaged_wake_uncropped, 'linestyle', 'none')
    plot(x, sgolayfilt(massaged_bottom_edge, 3, 21), 'color', 'black', 'linewidth', 2)
    colorbar()

    % Save to fit a line
    xs(block, :) = x;
    ys(block, :) = massaged_bottom_edge;
    
end


% Fit using all points
x_fit = reshape(xs.', 1, []);
y_fit = reshape(ys.', 1, []);


hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-3, 3])
xlabel('$x / D$', 'interpreter', 'latex')
ylabel('$z / D$', 'interpreter', 'latex')




%% Compute momentum thickness

u_inf = 7.5;
turbine = 1;
half_width_factor = 0.5;

figure('color', 'white')
hold on

% Store momentum thickness results
theta_raw = cell(1,3);
theta_massaged = cell(1,3);

for block = 1:3

    % Load coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y - turbine_positions(turbine);
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;

    % Crop inner side
    raw_wake(Y > 3) = nan;
    massaged_wake(Y > 3) = nan;

    % Get centerline velocities (core velocity U2(x))
    raw_center_velocity = mean_centerline_velocity(block).raw.u;         % 1 x Nx
    massaged_center_velocity = mean_centerline_velocity(block).massaged.u;

    % -----------------------------
    % Mixing-layer momentum thickness theta(x)
    % U1 = u_inf (outer), U2 = U_core(x) (inner/core)
    % theta(x) = ∫ ((U-U2)(U1-U) / (U1-U2)^2) dy
    % -----------------------------
    U1 = u_inf;

    % Preallocate
    theta_raw{block}      = nan(size(x));
    theta_massaged{block} = nan(size(x));

    % Column-wise integration over y
    for j = 1:numel(x)

        % --- RAW ---
        Ucol = raw_wake(:, j);
        U2   = raw_center_velocity(j);
        dU   = U1 - U2;

        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) >= 3 && isfinite(dU) && dU > 1e-6
            integrand = ((Ucol(good) - U2) .* (U1 - Ucol(good))) ./ (dU.^2);
            % Momentum thickness in rotor diameters D
            theta_raw{block}(j) = trapz(y(good), integrand);
        end

        % --- MASSAGED ---
        Ucol = massaged_wake(:, j);
        U2   = massaged_center_velocity(j);
        dU   = U1 - U2;

        good = isfinite(Ucol) & isfinite(y);
        if nnz(good) >= 3 && isfinite(dU) && dU > 1e-6
            integrand = ((Ucol(good) - U2) .* (U1 - Ucol(good))) ./ (dU.^2);
            % Momentum thickness in rotor diameters D
            theta_massaged{block}(j) = trapz(y(good), integrand);
        end
    end

    % Find half-width velocity
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
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);
    raw_bottom_edge = Y(raw_bottom_edge_index);
    massaged_bottom_edge = Y(massaged_bottom_edge_index);

    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y(1);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y(1);

    % Plot
    % contourf(X, Y, massaged_wake)
    % plot(x, sgolayfilt(massaged_bottom_edge, 3, 21), 'color', 'black', 'linewidth', 2)

    % Save to fit a line
    xs(block, :) = x;
    ys(block, :) = massaged_bottom_edge;

end

% Fit using all points
x_fit = reshape(xs.', 1, []);
y_fit = reshape(ys.', 1, []);

scatter(x_fit, y_fit, 20, 'filled')

hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-5, 6])
xlim([0, 51])
yline(0, 'linewidth', 1, 'color', 'black')







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
