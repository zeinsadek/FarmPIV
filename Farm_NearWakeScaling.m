% Attempt to show the near wake can be scaled
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
% SEGMENT TURBINES (NEAR WAKE, BLOCK 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turbine/Farm dimensions
D = 80;
Sy = 3;
u_inf = 7.5;

% Load data
block = 1;
X = data(block).X;
Y = data(block).Y;
x = X(1,:);
y = Y(:,1);

raw_u = data(block).u;
massaged_u = cleaned(block).u;

% T1 ~ edge turbine
% T4 ~ center turbine
turbine_positions = -9:3:3;

% Loop through and crop different wakes
for t = 1:length(turbine_positions)
    centered_y = y - turbine_positions(t);
    mask = (centered_y >= -Sy/2) & (centered_y <= Sy/2);

    % Save
    wakes(t).raw.u = raw_u(mask, :);
    wakes(t).massaged.u = massaged_u(mask, :);
    wakes(t).Y = Y(mask, :);
    wakes(t).X = X(mask, :);
end

clear t centered_y mask


clc; close all
% Compute center lines and half widths for 
half_width_factor = 0.5;
for t = 1:length(turbine_positions)

    % Raw + Massaged
    raw_wake = wakes(t).raw.u;
    massaged_wake = wakes(t).massaged.u;
    [r,~] = size(raw_wake);
    center_index = round(r/2);

    % Centerline of both datasets
    raw_center_velocity = raw_wake(center_index, :);
    massaged_center_velocity = massaged_wake(center_index, :);
    raw_center_half_velocity = half_width_factor * (raw_center_velocity + u_inf);
    massaged_center_half_velocity = half_width_factor * (massaged_center_velocity + u_inf);

    % Inside half-width
    raw_wake(raw_wake > raw_center_half_velocity) = nan;
    massaged_wake(massaged_wake > massaged_center_half_velocity) = nan;

    % Detect top and bottom half-widths
    raw_isnan = isnan(raw_wake);
    massaged_isnan = isnan(massaged_wake);



    %%% Top (bottom unflipped)
    Y_bottom = wakes(t).Y(center_index:end, 1);
    raw_isnan_bottom = raw_isnan(center_index:end, :);
    massaged_isnan_bottom = massaged_isnan(center_index:end, :);

    % Vertical difference to find edge
    raw_isnan_bottom_diff = abs(diff(raw_isnan_bottom, 1, 1));
    massaged_isnan_bottom_diff = abs(diff(massaged_isnan_bottom, 1, 1));
    [~, raw_bottom_edge_index] = max(raw_isnan_bottom_diff, [], 1);
    [~, massaged_bottom_edge_index] = max(massaged_isnan_bottom_diff, [], 1);
    raw_bottom_edge = Y_bottom(raw_bottom_edge_index);
    massaged_bottom_edge = Y_bottom(massaged_bottom_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_bottom_edge(raw_bottom_edge_index == 1) = Y_bottom(end);
    massaged_bottom_edge(massaged_bottom_edge_index == 1) = Y_bottom(end);


    %%% Bottom (top unflipped)
    Y_top = wakes(t).Y(1:center_index, 1);
    raw_isnan_top = raw_isnan(1:center_index, :);
    massaged_isnan_top = massaged_isnan(1:center_index, :);

    % Vertical difference to find edge
    raw_isnan_top_diff = abs(diff(raw_isnan_top, 1, 1));
    massaged_isnan_top_diff = abs(diff(massaged_isnan_top, 1, 1));
    [~, raw_top_edge_index] = max(raw_isnan_top_diff, [], 1);
    [~, massaged_top_edge_index] = max(massaged_isnan_top_diff, [], 1);
    raw_top_edge = Y_top(raw_top_edge_index);
    massaged_top_edge = Y_top(massaged_top_edge_index);
    
    % Fix values for when detected value goes out of bounds
    raw_top_edge(raw_top_edge_index == 1) = Y_top(1);
    massaged_top_edge(massaged_top_edge_index == 1) = Y_top(1);


    %%% Wake centerline
    wake_y = wakes(t).Y(:, 1);
    [~, raw_minimum_index] = min(raw_wake, [], 1);
    [~, massaged_minimum_index] = min(massaged_wake, [], 1);
    raw_center = wake_y(raw_minimum_index);
    massaged_center = wake_y(massaged_minimum_index);


    % Save center line velocity
    wakes(t).raw.center_line = raw_center_velocity;
    wakes(t).massaged.center_line = massaged_center_velocity;

    % Save half widths
    wakes(t).raw.bottom_edge = raw_bottom_edge;
    wakes(t).massaged.bottom_edge = massaged_bottom_edge;
    wakes(t).raw.top_edge = raw_top_edge;
    wakes(t).massaged.top_edge = massaged_top_edge;

    % Save wake center
    wakes(t).raw.centerline = raw_center;
    wakes(t).massaged.centerline = massaged_center;




    %%% Plotting
    figure('color', 'white')
    tiledlayout(1,2)
    sgtitle(sprintf('Turbine %1.0f', t))
    nexttile
    hold on
    contourf(wakes(t).X, wakes(t).Y, raw_wake, 100, 'linestyle', 'none')
    plot(x, raw_center, 'color', 'red', 'linewidth',1)
    plot(x, raw_bottom_edge, 'color', 'black', 'linewidth', 2)
    plot(x, raw_top_edge, 'color', 'black', 'linewidth', 2)
    hold off
    axis equal
    title('Raw')
    set(gca, 'YDir', 'reverse')
    xline([4,7])

    nexttile
    hold on
    contourf(wakes(t).X, wakes(t).Y, massaged_wake, 100, 'linestyle', 'none')
    plot(x, massaged_center, 'color', 'red', 'linewidth',1)
    plot(x, massaged_bottom_edge, 'color', 'black', 'linewidth', 2)
    plot(x, massaged_top_edge, 'color', 'black', 'linewidth', 2)
    hold off
    axis equal
    title('Massaged')
    set(gca, 'YDir', 'reverse')
    xline([4,7])
end

clear t r c center_half_velocity center_index center_half_velocity
clear block half_width_factor massaged_bottom_edge massaged_bottom_edge_index
clear massaged_center_half_velocity massaged_center_velocity massaged_isnan
clear massaged_isnan_bottom massaged_isnan_bottom_diff massaged_isnan_top
clear massaged_isnan_top_diff massaged_top_edge massaged_top_edge_index massaged_u
clear massaged_wake raw_bottom_edge raw_bottom_edge_index raw_center_half_velocity
clear raw_center_velocity raw_isnan raw_isnan_bottom raw_isnan_bottom_diff 
clear raw_isnan_top raw_isnan_top_diff raw_top_edge raw_top_edge_index 
clear raw_u raw_wake Y_bottom Y_top wake_y raw_minimum_index massaged_minimum_index
clear raw_center massaged_center


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHOW PROFILES COLLAPSE (NEAR WAKE, BLOCK 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we are taking the mean of the halfwidth on either side

% Fontsizes
labelFontSize = 12;

% Turbine colors
turbine_colors = slanCM('imola', 5);

% Near wake ranges from 1-10D
x_slice_positions = [1,3,5,6];

% Plot settings
data_type = 'raw';
marker_size = 10;
marker_alpha = 1;

clc; close all
figure('color', 'white')
t = tiledlayout(1,length(x_slice_positions));

% Loop through different x positions
for i = 1:length(x_slice_positions)

    x_slice_position = x_slice_positions(i);
    [~, index] = min(abs(x - x_slice_position));

    h(i) = nexttile;
    hold on
    title(sprintf('$x \\mathbin{/} D = %1.0f$', x_slice_position), 'interpreter', 'latex')

    % Reference gaussian
    eta = linspace(-4, 4, 400);                       % eta = y / r_{1/2}
    gauss_ref = exp(-log(2) * eta.^2);                % so gauss_ref(eta=1)=0.5
    plot(gauss_ref, eta, 'k-', 'LineWidth', 3, 'DisplayName', 'Gaussian ref')
   
    % Loop through different turbines
    for t = 1:length(turbine_positions)

        % Get turbine spanwise coordinate
        % turbine_y = wakes(t).Y(:,1) - turbine_positions(t);
        turbine_y = wakes(t).Y(:,1) - wakes(t).(data_type).centerline(index);

        % Get velocity profile
        u_profile = wakes(t).(data_type).u(:, index);
        u_deficit = u_inf - u_profile;
        u_deficit_normalized = u_deficit / (max(u_deficit, [], 'all', 'omitnan'));

        % Get half width (mean between the two sides, for now)
        top_half_width = wakes(t).(data_type).top_edge(index) - wakes(t).(data_type).centerline(index);
        bottom_haf_width = wakes(t).(data_type).bottom_edge(index) - wakes(t).(data_type).centerline(index);
        mean_half_width = mean([-top_half_width, bottom_haf_width], 'all', 'omitnan');

        % Plot
        scatter(u_deficit_normalized, turbine_y / mean_half_width, marker_size, 'filled', ...
                'MarkerFaceColor', turbine_colors(t,:), 'HandleVisibility', 'off', ...
                'MarkerFaceAlpha', marker_alpha)
        set(gca, 'YDir', 'Reverse')

        if i ~= 1
            ax = gca;
            ax.YAxis.Visible = 'off';
        end

        if i == 1
            ylabel('$z \mathbin{/} z_{1/2}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
        end
    end
end
hold off
linkaxes(h, 'xy')
xlim([-0.05, 1])
ylim([-4,4])

% Add turbine color legend
hold on
for t = 1:5
    scatter(nan, nan, 20, 'filled', 'MarkerFaceColor', turbine_colors(t,:), ...
            'DisplayName', sprintf('Turbine %1.0f', t))
end
hold off
leg = legend('interpreter', 'latex', 'orientation', 'horizontal', 'box', 'off');
leg.Layout.Tile = 'north';




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHOW PROFILES COLLAPSE (NEAR WAKE, BLOCK 1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here we are normalizing each half of the profile by the respective half
% width

% Fontsizes
labelFontSize = 12;

% Turbine colors
turbine_colors = slanCM('imola', 5);

% Near wake ranges from 1-10D
x_slice_positions = [1,3,5,6];

% Plot settings
data_type = 'raw';
marker_size = 10;
marker_alpha = 1;

clc; close all
figure('color', 'white')
t = tiledlayout(1,length(x_slice_positions));

% Loop through different x positions
for i = 1:length(x_slice_positions)

    x_slice_position = x_slice_positions(i);
    [~, index] = min(abs(x - x_slice_position));

    h(i) = nexttile;
    hold on
    title(sprintf('$x \\mathbin{/} D = %1.0f$', x_slice_position), 'interpreter', 'latex')

    % Reference gaussian
    eta = linspace(-4, 4, 400);                       % eta = y / r_{1/2}
    gauss_ref = exp(-log(2) * eta.^2);                % so gauss_ref(eta=1)=0.5
    plot(gauss_ref, eta, 'k-', 'LineWidth', 3, 'DisplayName', 'Gaussian ref')
   
    % Loop through different turbines
    for t = 1:length(turbine_positions)

        % Get turbine spanwise coordinate
        turbine_y = wakes(t).Y(:,1) - wakes(t).(data_type).centerline(index);

        % Get velocity profile
        u_profile = wakes(t).(data_type).u(:, index);
        u_deficit = u_inf - u_profile;
        u_deficit_normalized = u_deficit / (max(u_deficit, [], 'all', 'omitnan'));

        % Get half width (mean between the two sides, for now)
        top_half_width = abs(wakes(t).(data_type).top_edge(index) - wakes(t).(data_type).centerline(index));
        bottom_haf_width = abs(wakes(t).(data_type).bottom_edge(index) - wakes(t).(data_type).centerline(index));

        % Normalize by half-widths
        center_index = round(length(turbine_y)/2);
        normalized_turbine_y = turbine_y;
        normalized_turbine_y(1:center_index) = turbine_y(1:center_index) / bottom_haf_width;
        normalized_turbine_y(center_index:end) = turbine_y(center_index:end) / top_half_width;


        % Plot
        scatter(u_deficit_normalized, normalized_turbine_y, marker_size, 'filled', ...
                'MarkerFaceColor', turbine_colors(t,:), 'HandleVisibility', 'off', ...
                'MarkerFaceAlpha', marker_alpha)
        set(gca, 'YDir', 'Reverse')

        if i ~= 1
            ax = gca;
            ax.YAxis.Visible = 'off';
        end

        if i == 1
            ylabel('$z \mathbin{/} z_{1/2}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
        end
    end
end
hold off
linkaxes(h, 'xy')
xlim([-0.05, 1])
ylim([-4,4])

% Add turbine color legend
hold on
for t = 1:5
    scatter(nan, nan, 20, 'filled', 'MarkerFaceColor', turbine_colors(t,:), ...
            'DisplayName', sprintf('Turbine %1.0f', t))
end
hold off
leg = legend('interpreter', 'latex', 'orientation', 'horizontal', 'box', 'off');
leg.Layout.Tile = 'north';



%% Compare edge turbine profiles to the center turbine
% Try modeling as a skewnormal profile



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTERLINE VELOCITIES (COMPARE DIFFERENT DEFINITIONS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp = nan(5, 195);

clear h

clc; close all
figure('color', 'white')
title('Near Wake')
hold on
for t = 1:5
    raw_wake = wakes(t).raw.u;
    massaged_wake = wakes(t).massaged.u;
    wake_X = wakes(t).X;
    wake_Y = wakes(t).Y - turbine_positions(t);
    wake_x = wake_X(1,:);
    wake_y = wake_Y(:,1);

    % Get centerline by drawing a line straight back
    [r,~] = size(raw_wake);
    center_index = round(r/2);

    % Centerline straight back
    raw_center_velocity = raw_wake(center_index, :);
    massaged_center_velocity = massaged_wake(center_index, :);
    
    % Normalize
    raw_center_velocity_deficit = u_inf - raw_center_velocity;
    massaged_center_velocity_deficit = u_inf - massaged_center_velocity;

    norm_raw_center_velocity_deficit = raw_center_velocity_deficit / u_inf;
    norm_massaged_center_velocity_deficit = massaged_center_velocity_deficit / u_inf;



    % Centerline by following wake
    raw_wake_following_center_velocity = min(raw_wake, [], 1, 'omitnan');
    massaged_wake_following_center_velocity = min(massaged_wake, [], 1, 'omitnan');

    % Normalize
    raw_wake_following_center_velocity_deficit = u_inf - raw_wake_following_center_velocity;
    massaged_wake_following_center_velocity_deficit = u_inf - massaged_wake_following_center_velocity;

    norm_raw_wake_following_center_velocity_deficit = raw_wake_following_center_velocity_deficit / u_inf;
    norm_massaged_wake_following_center_velocity_deficit = massaged_wake_following_center_velocity_deficit / u_inf;


    % Plot all centerlines together
    plot(wake_x, norm_massaged_center_velocity_deficit, ...
         'color', turbine_colors(t, :), 'linewidth', 2, ...
         'DisplayName', sprintf('Turbine %1.0f', t))

    % Collect all profiles to average over turbines
    tmp(t, :) = norm_massaged_center_velocity_deficit;


    % % Plot to compare
    % figure('color', 'white')
    % tiledlayout(1,2)
    % sgtitle(sprintf('Turbine %1.0f Center line velocities', t))
    % 
    % % Raw
    % h(1) = nexttile;
    % hold on
    % plot(wake_x, norm_raw_center_velocity_deficit, 'linestyle', '-')
    % plot(wake_x, norm_raw_wake_following_center_velocity_deficit, 'linestyle', '--')
    % hold off
    % title('Raw')
    % 
    % % Massaged
    % h(1) = nexttile;
    % hold on
    % plot(wake_x, norm_massaged_center_velocity_deficit, 'linestyle', '-')
    % plot(wake_x, norm_massaged_wake_following_center_velocity_deficit, 'linestyle', '--')
    % hold off
    % title('Massaged')
    % linkaxes(h, 'xy')
    % ylim([0, 1])

end

% Plot average over all turbines
mean_norm_centerline_deficit = mean(tmp, 1, 'omitnan');
plot(wake_x, mean_norm_centerline_deficit, 'color', 'black', 'linewidth', 3, 'DisplayName', 'Mean')
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')

% JENSEN: Out-of-the-box
Ct = 0.644;
k = 0.07;
A0 = 1 - sqrt(1 - Ct);
Jensen_amplitude = (1 + (2 * k * wake_x)).^(-2);
% plot(wake_x, Jensen_amplitude, 'linewidth', 2, 'color', 'blue', 'DisplayName', 'Jensen');


% JENSEN: Optimized
% xD = wake_x;                       % this is x/D in your plotting
% 
% xD_min = 1.;   % exclude near-wake (tweak: 1D-3D are common cutoffs)
% xD_max = 10.0;  % match your axis limit or choose smaller
% 
% mask = isfinite(xD) & isfinite(mean_norm_centerline_deficit) & (xD >= xD_min) & (xD <= xD_max);

xFit = wake_x;
yFit = mean_norm_centerline_deficit;

% Jensen normalized deficit shape
jensen_norm = @(k, x) (1 + 2*k.*x).^(-2);
obj = @(k) sum((yFit - jensen_norm(k, xFit)).^2);
k0 = 0.07;

% fminsearch is unconstrained; enforce k>0 via reparameterization
% k = exp(p) ensures positivity
obj_p = @(p) sum( (yFit - jensen_norm(exp(p), xFit)).^2 );
p0 = log(k0);
p_hat = fminsearch(obj_p, p0);
k_hat = exp(p_hat);

% Best-fit curve on full x range for plotting
yHat = jensen_norm(k_hat, wake_x);

% Fit diagnostics
SSE = sum((yFit - jensen_norm(k_hat, xFit)).^2);
SST = sum((yFit - mean(yFit)).^2);
R2  = 1 - SSE/SST;

% Print results to command window
% fprintf('Jensen best-fit k = %.5f (fit range x/D = %.2f to %.2f)\n', k_hat, xD_min, xD_max);
% fprintf('SSE = %.4g, R^2 = %.4f, N = %d\n', SSE, R2, numel(xFit));

% Optional: annotate on the plot
txt = sprintf('Best-fit Jensen: k = %.3f, R^2 = %.3f', k_hat, R2);
% text(0.55*xD_max, 0.15, txt, 'Color', 'k', 'FontWeight', 'bold');


% Plot
plot(wake_x, yHat, 'color', 'red', 'linestyle', '--',  'LineWidth', 2, 'DisplayName', sprintf('Jensen Optimized\n$k = %1.3f$', k_hat));
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')








% BPA14 fit (DeltaU/Uinf)_max: optimize k*, epsilon, Ct WITHOUT clipping

xD = wake_x;        % x/D
y  = mean_norm_centerline_deficit;   % (DeltaU/Uinf)_max  (dimensionless)

% --- fit window (strongly suggest excluding x/D < ~2 for BPA) ---
xD_min = 1.0;
xD_max = 10.0;

mask = isfinite(xD) & isfinite(y) & (xD >= xD_min) & (xD <= xD_max);
xFit = xD(mask);
yFit = y(mask);

% ----- transforms: k*>0, eps>0, 0<Ct<1 -----
invlogit = @(z) 1 ./ (1 + exp(-z));
logit    = @(c) log(c ./ (1-c));

% Model (no clipping)
C_model = @(kstar, eps, Ct, x) 1 - sqrt( 1 - Ct ./ (8*(kstar.*x + eps).^2) );

% Objective with validity penalty (no NaN clipping)
obj = @(p) obj_bpa_valid(p, xFit, yFit, C_model, invlogit);

% Initial guesses
k0   = 0.05;
eps0 = 0.20;
Ct0  = 0.644;

p0 = [log(k0), log(eps0), logit(Ct0)];

% Fit
opts = optimset('Display','iter');   % show progress so you can see if it's moving
p_hat = fminsearch(obj, p0, opts);
k_hat  = exp(p_hat(1));
eps_hat = exp(p_hat(2));
Ct_hat  = invlogit(p_hat(3));

% Diagnostics
yHat_fit = C_model(k_hat, eps_hat, Ct_hat, xFit);
SSE = sum((yFit - yHat_fit).^2, 'omitnan');
SST = sum((yFit - mean(yFit,'omitnan')).^2, 'omitnan');
R2  = 1 - SSE/SST;

fprintf('BPA14 best-fit: k* = %.5f, eps = %.5f, Ct = %.5f, R^2 = %.4f\n', k_hat, eps_hat, Ct_hat, R2);

% Plot
yHat_full = C_model(k_hat, eps_hat, Ct_hat, xD);
bastankhah_label = sprintf('Bastankhah Optimized\n$C_t = %1.3f$\n$k = %1.3f$\n$\\varepsilon = %1.3f$', Ct_hat, k_hat, eps_hat);

plot(xD, yHat_full, 'color', 'green', 'linestyle', '--', 'LineWidth', 2, 'DisplayName', bastankhah_label);
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% Ishihara & Qian-style fit:  A(x̃) = (a + b x̃ + c (1+x̃)^(-2))^(-2)

xT = wake_x;        % x̃ = x/D
y  = mean_norm_centerline_deficit;   % (ΔU/U∞)_max

% --- fit window ---
xT_min = 1.0;
xT_max = 10.0;

mask = isfinite(xT) & isfinite(y) & (xT >= xT_min) & (xT <= xT_max);
xFit = xT(mask);
yFit = y(mask);

% Model
A_model = @(p, x) (p(1) + p(2).*x + p(3).*(1+x).^(-2)).^(-2);   % p=[a b c]
g_model = @(p, x) (p(1) + p(2).*x + p(3).*(1+x).^(-2));        % denominator g(x)

% Objective: SSE + penalty if g<=0 anywhere
obj = @(p) obj_IQ(p, xFit, yFit, A_model, g_model);

% Initial guess (reasonable-ish)
% If y ~ 0.6 at x=1, then g ~ 1/sqrt(y) ~ 1.29. Use that as a ballpark.
p0 = [1.2, 0.05, 0.5];   % [a b c]

opts = optimset('Display','iter');
p_hat = fminsearch(obj, p0, opts);

a_hat = p_hat(1);
b_hat = p_hat(2);
c_hat = p_hat(3);


% Diagnostics
yHat_fit = A_model(p_hat, xFit);
SSE = sum((yFit - yHat_fit).^2, 'omitnan');
SST = sum((yFit - mean(yFit,'omitnan')).^2, 'omitnan');
R2  = 1 - SSE/SST;

fprintf('Ishihara-Qian fit: a=%.5f, b=%.5f, c=%.5f, R^2=%.4f\n', a_hat, b_hat, c_hat, R2);

% Evaluate + plot
yHat_full = A_model(p_hat, xT);

ishihara_label = sprintf('Ishihara Optimized\n$a = %1.3f$\n$b = %1.3f$\n$c = %1.3f$', a_hat, b_hat, c_hat);
plot(xT, yHat_full, 'color', 'magenta', 'linestyle', '--', 'LineWidth', 2, 'DisplayName', ishihara_label);













hold off
xlim([0, 10])
ylim([0, 1.1])
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$\frac{\Delta u}{u_{\infty}}$', 'interpreter', 'latex')
legend('interpreter', 'latex', 'box', 'off', 'location', 'northeast')










%% Functions

% ---- local objective with validity penalty ----
function J = obj_bpa_valid(p, x, y, C_model, invlogit)

    k   = exp(p(1));
    eps = exp(p(2));
    Ct  = invlogit(p(3));

    s = k.*x + eps;                         % sigma/d0 = k*x + eps
    arg = 1 - Ct ./ (8*s.^2);               % inside sqrt

    % Penalty if arg goes negative anywhere in fit window
    minArg = min(arg);
    if minArg <= 0
        % Quadratic penalty grows as you violate constraint
        J = 1e6 + 1e8*(minArg)^2;
        return
    end

    yModel = C_model(k, eps, Ct, x);

    if any(~isfinite(yModel)) || any(~isreal(yModel))
        J = 1e9;
        return
    end

    r = y - yModel;
    J = sum(r.^2, 'omitnan');

end





% --- local objective ---
function J = obj_IQ(p, x, y, A_model, g_model)

    g = g_model(p, x);

    % Hard penalty if denominator nonpositive (prevents blow-ups/complex)
    minG = min(g);
    if minG <= 0
        J = 1e6 + 1e8*(minG)^2;
        return
    end

    yModel = A_model(p, x);
    if any(~isfinite(yModel)) || any(~isreal(yModel))
        J = 1e9;
        return
    end

    r = y - yModel;
    J = sum(r.^2, 'omitnan');

    % Optional: mild regularization to avoid absurd parameters
    % J = J + 1e-4*sum(p.^2);

end



