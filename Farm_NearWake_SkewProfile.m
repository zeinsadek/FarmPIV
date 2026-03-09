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




    % %%% Plotting
    % figure('color', 'white')
    % tiledlayout(1,2)
    % sgtitle(sprintf('Turbine %1.0f', t))
    % nexttile
    % hold on
    % contourf(wakes(t).X, wakes(t).Y, raw_wake, 100, 'linestyle', 'none')
    % plot(x, raw_center, 'color', 'red', 'linewidth',1)
    % plot(x, raw_bottom_edge, 'color', 'black', 'linewidth', 2)
    % plot(x, raw_top_edge, 'color', 'black', 'linewidth', 2)
    % hold off
    % axis equal
    % title('Raw')
    % set(gca, 'YDir', 'reverse')
    % xline([4,7])
    % 
    % nexttile
    % hold on
    % contourf(wakes(t).X, wakes(t).Y, massaged_wake, 100, 'linestyle', 'none')
    % plot(x, massaged_center, 'color', 'red', 'linewidth',1)
    % plot(x, massaged_bottom_edge, 'color', 'black', 'linewidth', 2)
    % plot(x, massaged_top_edge, 'color', 'black', 'linewidth', 2)
    % hold off
    % axis equal
    % title('Massaged')
    % set(gca, 'YDir', 'reverse')
    % xline([4,7])
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
    % for t = 1:length(turbine_positions)
    for t = 1:3:4

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
for t = 1:3:4
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
    % for t = 1:length(turbine_positions)
    for t = 1:3:4

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


% Fontsizes
labelFontSize = 12;

% Turbine colors
turbine_colors = slanCM('imola', 5);

% Near wake ranges from 1-10D
x_slice_positions = [1,3,5,6,6.8,8];

% Plot settings
data_type = 'raw';
marker_size = 10;
marker_alpha = 1;

clc; close all
figure('color', 'white')
tiledlayout(1,length(x_slice_positions));

% Loop through different x positions
for i = 1:length(x_slice_positions)

    x_slice_position = x_slice_positions(i);
    [~, index] = min(abs(x - x_slice_position));

    h(i) = nexttile;
    hold on
    title(sprintf('$x \\mathbin{/} D = %1.1f$', x_slice_position), 'interpreter', 'latex')

    for t = 1:3:4

        % Get velocity profile
        u_profile = wakes(t).(data_type).u(:, index);
        u_deficit = u_inf - u_profile;
        u_deficit = u_deficit / u_inf;
        % u_deficit_normalized = u_deficit / u_inf;
        u_deficit_normalized = u_deficit / (max(u_deficit, [], 'all', 'omitnan'));
        
        % Get turbine spanwise coordinate
        turbine_y = wakes(t).Y(:,1) - turbine_positions(t);

        % Plot
        scatter(u_deficit_normalized, turbine_y, marker_size, 'filled', ...
                'MarkerFaceColor', turbine_colors(t,:), 'HandleVisibility', 'off', ...
                'MarkerFaceAlpha', marker_alpha)
        set(gca, 'YDir', 'Reverse')

        if i ~= 1
            ax = gca;
            ax.YAxis.Visible = 'off';
        end

        if i == 1
            ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
        end
    end
end

hold off
linkaxes(h, 'xy')
xlim([-0.01, 1])
ylim([-1.5,1.5])

% Add turbine color legend
hold on
for t = 1:3:4
    scatter(nan, nan, 20, 'filled', 'MarkerFaceColor', turbine_colors(t,:), ...
            'DisplayName', sprintf('Turbine %1.0f', t))
end
hold off
leg = legend('interpreter', 'latex', 'orientation', 'horizontal', 'box', 'off');
leg.Layout.Tile = 'north';





%% Compare edge turbine profiles to the center turbine
% Fit symmetric Gaussian and skew-normal to deficit vs y

labelFontSize = 12;
turbine_colors = slanCM('imola', 5);

x_slice_positions = [1,3,5,7,9];

data_type = 'raw';
marker_size = 10;
marker_alpha = 1;

clc; close all
figure('color', 'white')
tl = tiledlayout(1, length(x_slice_positions));

% --- Models
gaussModel = @(p,y) p(1) * exp(-(y - p(2)).^2 ./ (2*p(3)^2)) + p(4);                 % [A,mu,sigma,B]
skewModel  = @(p,y) p(1) .* (2./p(3)) .* normpdf((y - p(2))./p(3)) .* ...
                           normcdf(p(4).*(y - p(2))./p(3)) + p(5);                  % [A,xi,omega,alpha,B]

% Fit options
opts = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',2e4);

% Store results if you want
fitResults = struct();

axHandles = gobjects(1, length(x_slice_positions));

for i = 1:length(x_slice_positions)

    x_slice_position = x_slice_positions(i);
    [~, index] = min(abs(x - x_slice_position));

    axHandles(i) = nexttile;
    hold on
    title(sprintf('$x \\mathbin{/} D = %1.1f$', x_slice_position), 'interpreter', 'latex')

    for tt = 1:3:4  % (rename from t -> tt)

        % Get velocity profile (as function of y)
        u_profile = wakes(tt).(data_type).u(:, index);
        u_deficit = (u_inf - u_profile) ./ u_inf;

        % Normalize by its max (what you're currently doing)
        umax = max(u_deficit, [], 'all', 'omitnan');
        u_deficit_normalized = u_deficit ./ umax;

        % Spanwise coordinate (your variable name says y, but label says z/D — keep consistent with your convention)
        turbine_y = wakes(tt).Y(:,1) - turbine_positions(tt);

        % --- Clean data for fitting
        ydat = turbine_y(:);
        xdat = u_deficit_normalized(:);

        good = isfinite(ydat) & isfinite(xdat);
        ydat = double(ydat(good));
        xdat = double(xdat(good));

        % Optional: restrict to a window to avoid far-field noise
        % keep = abs(ydat) <= 1.5;
        % ydat = ydat(keep); xdat = xdat(keep);

        % Plot raw points
        scatter(xdat, ydat, marker_size, 'filled', ...
            'MarkerFaceColor', turbine_colors(tt,:), 'HandleVisibility', 'off', ...
            'MarkerFaceAlpha', marker_alpha)

        % --- Initial guesses
        [xpk, idxpk] = max(xdat);
        mu0    = ydat(idxpk);
        sig0   = 0.25;              % tweak if needed
        A0     = max(xpk, 0.8);      % near 1 usually
        B0     = max(min(xdat), 0);  % baseline near 0
        % Gaussian init
        p0g = [A0, mu0, sig0, B0];

        % Skew init: start with small skew
        alpha0 = 0.0;
        omega0 = sig0;
        p0s = [A0, mu0, omega0, alpha0, B0];

        % --- Bounds (keep things sane)
        % Gaussian bounds
        lbg = [0,   min(ydat),  0.02, -0.2];
        ubg = [2.0, max(ydat),  2.00,  0.5];

        % Skew bounds
        lbs = [0,   min(ydat),  0.02, -20,  -0.2];
        ubs = [2.0, max(ydat),  2.00,  20,   0.5];

        % --- Fit
        try
            pg = lsqcurvefit(gaussModel, p0g, ydat, xdat, lbg, ubg, opts);
            xsse_g = sum((xdat - gaussModel(pg, ydat)).^2);

            ps = lsqcurvefit(skewModel,  p0s, ydat, xdat, lbs, ubs, opts);
            xsse_s = sum((xdat - skewModel(ps, ydat)).^2);
        catch ME
            warning('Fit failed at x/D=%g, turbine=%d: %s', x_slice_position, tt, ME.message);
            continue
        end

        % --- Overlay fitted curves
        yline = linspace(min(ydat), max(ydat), 400);

        xg = gaussModel(pg, yline);
        xs = skewModel(ps, yline);

        plot(xg, yline, '-',  'LineWidth', 2, 'Color', turbine_colors(tt,:), 'HandleVisibility','off')
        plot(xs, yline, '--', 'LineWidth', 2, 'Color', turbine_colors(tt,:), 'HandleVisibility','off')

        % --- Save results
        fitResults(i).xD = x_slice_position;
        fitResults(i).turbine(tt).gauss.p = pg;
        fitResults(i).turbine(tt).skew.p  = ps;
        fitResults(i).turbine(tt).sse_gauss = xsse_g;
        fitResults(i).turbine(tt).sse_skew  = xsse_s;

        % Optional: quick "which is better" metric (AIC-ish)
        n = numel(xdat);
        kg = 4; ks = 5;
        AICg = n*log(xsse_g/n) + 2*kg;
        AICs = n*log(xsse_s/n) + 2*ks;
        fitResults(i).turbine(tt).AIC_gauss = AICg;
        fitResults(i).turbine(tt).AIC_skew  = AICs;

    end

    set(gca, 'YDir', 'Reverse')

    if i ~= 1
        ax = gca;
        ax.YAxis.Visible = 'off';
    else
        ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
    end
end

linkaxes(axHandles, 'xy')
xlim([-0.01, 1])
ylim([-1.5, 1.5])

% Legend for turbines + line styles
hold on
for tt = 1:3:4
    scatter(nan, nan, 20, 'filled', 'MarkerFaceColor', turbine_colors(tt,:), ...
        'DisplayName', sprintf('Turbine %1.0f', tt))
end
plot(nan, nan, 'k-',  'LineWidth', 2, 'DisplayName', 'Gaussian fit')
plot(nan, nan, 'k--', 'LineWidth', 2, 'DisplayName', 'Skew-normal fit')
hold off

leg = legend('interpreter', 'latex', 'orientation', 'horizontal', 'box', 'off');
leg.Layout.Tile = 'north';






%% Batch fit wake profiles at each x/D and plot fit parameters vs x/D
% Fits: (1) symmetric Gaussian, (2) skew-normal
% Outputs summary plots vs x/D for turbines 1 and 4

% ---------------- User settings ----------------
data_type = 'massaged';
% xD_list   = [1,2,3,4,5,6,7,8,9,10];     % downstream locations (x/D)
xD_list   = 1:0.5:10;     % downstream locations (x/D)
turbines  = [1,4];               % center vs edge (adjust if needed)
y_window  = 1.5;                 % keep |y|<=y_window (set [] to disable)
minPts    = 25;                  % minimum points needed to attempt fit

% ---------------- Models ----------------
gaussModel = @(p,y) p(1) * exp(-(y - p(2)).^2 ./ (2*p(3)^2)) + p(4);  % [A,mu,sigma,B]

skewModel  = @(p,y) p(1) .* (2./p(3)) .* normpdf((y - p(2))./p(3)) .* ...
                           normcdf(p(4).*(y - p(2))./p(3)) + p(5);   % [A,xi,omega,alpha,B]

opts = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',2e4);

% ---------------- Storage arrays ----------------
nx = numel(xD_list);
nt = numel(turbines);

% Gaussian params
A_g    = nan(nx,nt);
mu_g   = nan(nx,nt);
sig_g  = nan(nx,nt);
B_g    = nan(nx,nt);
SSE_g  = nan(nx,nt);
AIC_g  = nan(nx,nt);

% Skew params
A_s     = nan(nx,nt);
xi_s    = nan(nx,nt);
om_s    = nan(nx,nt);
alpha_s = nan(nx,nt);
B_s     = nan(nx,nt);
SSE_s   = nan(nx,nt);
AIC_s   = nan(nx,nt);

% Convenience
idx_x = nan(nx,1);

% ---------------- Main fitting loop (NO plotting) ----------------
for i = 1:nx
    xD = xD_list(i);
    [~, idx_x(i)] = min(abs(x - xD));
    col = idx_x(i);

    for j = 1:nt
        tt = turbines(j);

        % Profile vs y
        u_profile = wakes(tt).(data_type).u(:, col);
        u_deficit = (u_inf - u_profile) ./ u_inf;

        umax = max(u_deficit, [], 'omitnan');
        if ~isfinite(umax) || umax <= 0
            continue
        end
        xdat = u_deficit ./ umax;

        ydat = wakes(tt).Y(:,1) - turbine_positions(tt);

        % Clean
        good = isfinite(xdat) & isfinite(ydat);
        xdat = double(xdat(good));
        ydat = double(ydat(good));

        % Optional windowing
        if ~isempty(y_window)
            keep = abs(ydat) <= y_window;
            xdat = double(xdat(keep));
            ydat = double(ydat(keep));
        end

        if numel(xdat) < minPts
            continue
        end

        % Initial guesses (peak-based)
        [xpk, kpk] = max(xdat);
        mu0   = ydat(kpk);
        sig0  = 0.25;                 % tweak if needed
        A0    = max(xpk, 0.8);
        B0    = max(min(xdat), 0);

        % Initial params
        p0g = [A0, mu0, sig0, B0];
        p0s = [A0, mu0, sig0, 0.0, B0];

        % Bounds
        lbg = [0,   min(ydat), 0.02, -0.2];
        ubg = [2.0, max(ydat), 2.00,  0.5];

        lbs = [0,   min(ydat), 0.02, -20, -0.2];
        ubs = [2.0, max(ydat), 2.00,  20,  0.5];

        % Fit
        try
            pg = lsqcurvefit(gaussModel, p0g, ydat, xdat, lbg, ubg, opts);
            resg = xdat - gaussModel(pg, ydat);
            sseg = sum(resg.^2);

            ps = lsqcurvefit(skewModel,  p0s, ydat, xdat, lbs, ubs, opts);
            ress = xdat - skewModel(ps, ydat);
            sses = sum(ress.^2);

        catch
            continue
        end

        % Save Gaussian
        A_g(i,j)   = pg(1);
        mu_g(i,j)  = pg(2);
        sig_g(i,j) = pg(3);
        B_g(i,j)   = pg(4);
        SSE_g(i,j) = sseg;

        % Save Skew
        A_s(i,j)     = ps(1);
        xi_s(i,j)    = ps(2);
        om_s(i,j)    = ps(3);
        alpha_s(i,j) = ps(4);
        B_s(i,j)     = ps(5);
        SSE_s(i,j)   = sses;

        % AIC (AIC-ish; good for relative comparison at same x/D)
        n  = numel(xdat);
        kg = 4;
        ks = 5;

        % Guard against log(0)
        sseg = max(sseg, eps);
        sses = max(sses, eps);

        AIC_g(i,j) = n*log(sseg/n) + 2*kg;
        AIC_s(i,j) = n*log(sses/n) + 2*ks;
    end
end

dAIC = AIC_g - AIC_s;   % positive => skew is better

% ---------------- Summary plots vs x/D ----------------

clc; close all
% 1) Skewness parameter alpha (key story)
figure('color','white'); hold on
plot(xD_list, alpha_s(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', turbines(1)), ...
    'color', turbine_colors(turbines(1), :), 'markerfacecolor', turbine_colors(turbines(1), :))
plot(xD_list, alpha_s(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', turbines(2)), ...
    'color', turbine_colors(turbines(2), :), 'markerfacecolor', turbine_colors(turbines(2), :))

xlabel('$x/D$','Interpreter','latex')
ylabel("Skew parameter $\alpha$",'Interpreter','latex')
legend('Interpreter', 'latex', 'Location', 'northwest', 'box', 'off', 'fontsize', 10); box off

% 2) ΔAIC: does skew-normal beat Gaussian?
% figure('color','white'); hold on
% plot(xD_list, dAIC(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', turbines(1)))
% plot(xD_list, dAIC(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', turbines(2)))
% xlabel('$x/D$','Interpreter','latex')
% ylabel('$\Delta$AIC = AIC$_G$ - AIC$_{skew}$','Interpreter','latex')
% legend('Interpreter','latex','Location','best'); box off

% 3) Center location vs x/D (mu / xi)
% figure('color','white'); hold on
% plot(xD_list, mu_g(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Gaussian $\mu$, T%d', turbines(1)))
% plot(xD_list, mu_g(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Gaussian $\mu$, T%d', turbines(2)))
% plot(xD_list, xi_s(:,1), '--o', 'LineWidth', 2, 'DisplayName', sprintf('Skew $\xi$, T%d', turbines(1)))
% plot(xD_list, xi_s(:,2), '--o', 'LineWidth', 2, 'DisplayName', sprintf('Skew $\xi$, T%d', turbines(2)))
% xlabel('$x/D$','Interpreter','latex')
% ylabel('Center location (in $D$ units)','Interpreter','latex')
% legend('Interpreter','latex','Location','best'); box off

% 4) Width vs x/D (sigma / omega)
% figure('color','white'); hold on
% plot(xD_list, sig_g(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Gaussian $\sigma$, T%d', turbines(1)))
% plot(xD_list, sig_g(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Gaussian $\sigma$, T%d', turbines(2)))
% plot(xD_list, om_s(:,1), '--o', 'LineWidth', 2, 'DisplayName', sprintf('Skew $\omega$, T%d', turbines(1)))
% plot(xD_list, om_s(:,2), '--o', 'LineWidth', 2, 'DisplayName', sprintf('Skew $\omega$, T%d', turbines(2)))
% xlabel('$x/D$','Interpreter','latex')
% ylabel('Width parameter','Interpreter','latex')
% legend('Interpreter','latex','Location','best'); box off

% (Optional) 5) SSE comparison (sometimes less interpretable than AIC)
figure('color','white'); hold on
plot(xD_list, SSE_g(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d: Gaussian', turbines(1)), ...
    'color', turbine_colors(turbines(1), :), 'markerfacecolor', turbine_colors(turbines(1), :))
plot(xD_list, SSE_s(:,1), '--o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d: Skew', turbines(1)), ...
    'color', turbine_colors(turbines(1), :), 'markerfacecolor', turbine_colors(turbines(1), :))

plot(nan, nan, 'color', 'white', 'displayname', ' ')

plot(xD_list, SSE_g(:,2), '-s', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d: Gaussian', turbines(2)), ...
    'color', turbine_colors(turbines(2), :), 'markerfacecolor', turbine_colors(turbines(2), :))
plot(xD_list, SSE_s(:,2), '--s', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d: Skew', turbines(2)), ...
    'color', turbine_colors(turbines(2), :), 'markerfacecolor', turbine_colors(turbines(2), :))
xlabel('$x/D$','Interpreter','latex')
ylabel('SSE','Interpreter','latex')
legend('Interpreter', 'latex', 'Location', 'northwest', 'box', 'off', 'fontsize', 10); box off














%% Model-free asymmetry metrics vs x/D (halfway-between-turbines cropping)

data_type = 'massaged';

% Turbines: (as you clarified)
tEdge   = 1;
tCenter = 4;
turbines = [tEdge, tCenter];

% Downstream sampling
xD_list = 1:0.5:10;  % or your custom [1,3,5,6,6.8,8] etc.

% Cropping bounds: half spacing = 1.5D if spacing is 3D
halfSpacing = 1.5;
bounds_mode = struct('type','window','zmin',-halfSpacing,'zmax',+halfSpacing);

% Choose deficit definition:
% A) shape-only: d = ΔU/ΔU_max  (recommended if comparing skew/shape)
% B) physical deficit: d = ΔU/U_inf
useMaxNorm = true;

% Optional core threshold (works best with max-normalized)
% Set [] to disable. Try 0.05–0.15.
coreThresh = 0.5;

% Where to split the halves
z0_mode = 'peak';   % 'peak' is robust (center at max deficit)

% Allocate storage
nx = numel(xD_list);
nt = numel(turbines);

I   = nan(nx,nt);   % signed imbalance
Rr  = nan(nx,nt);   % ratio >=1
dzc = nan(nx,nt);   % centroid shift relative to z0
z0  = nan(nx,nt);   % split location
At  = nan(nx,nt);   % total area (optional)

for i = 1:nx
    xD = xD_list(i);
    [~, col] = min(abs(x - xD));

    for j = 1:nt
        tt = turbines(j);

        % Pull profile at this x/D
        u_profile = wakes(tt).(data_type).u(:, col);

        % Spanwise coordinate in D units (you call it z/D)
        zD = wakes(tt).Y(:,1) - turbine_positions(tt);

        % Deficit
        d = (u_inf - u_profile) ./ u_inf;     % ΔU/U_inf

        if useMaxNorm
            dmax = max(d, [], 'omitnan');
            if ~isfinite(dmax) || dmax <= 0
                continue
            end
            d = d ./ dmax;                   % ΔU/ΔU_max
            thr = coreThresh;                % threshold makes sense here
        else
            thr = [];                        % threshold optional; often omit here
        end

        % Compute metrics inside +/- halfSpacing
        S = wake_asymmetry_metrics(zD, d, z0_mode, bounds_mode, ...
                                  'Threshold', thr, 'MinPts', 10);

        if ~S.valid
            continue
        end

        I(i,j)   = S.I;
        Rr(i,j)  = S.R;
        dzc(i,j) = S.dzc;
        z0(i,j)  = S.z0;
        At(i,j)  = S.Atot;
    end
end

% Plot: signed imbalance I vs x/D
figure('color','white'); hold on
plot(xD_list, I(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', tEdge), ...
    'color', turbine_colors(turbines(1), :), 'markerfacecolor', turbine_colors(turbines(1), :))

plot(xD_list, I(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Turbine %d', tCenter), ...
    'color', turbine_colors(turbines(2), :), 'markerfacecolor', turbine_colors(turbines(2), :))

set(gca, 'YDir', 'reverse')
xlabel('$x/D$','Interpreter','latex')
ylabel('$I = (A_R-A_L)/(A_R+A_L)$','Interpreter','latex')
legend('Interpreter','latex','Location','northwest', 'box', 'off'); box off
set(gca, 'YDir', 'reverse')

% Shade wake-merging region (edit limits as needed)
ylim([-0.4, 0.4])
mergeL = 5.5; mergeR = 7.0;
yl = ylim;
patch([mergeL mergeR mergeR mergeL],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 0 0], 'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
uistack(findobj(gca,'Type','line'),'top')

%% Plot: centroid shift dzc vs x/D
figure('color','white'); hold on
plot(xD_list, dzc(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Edge (T%d)', tEdge));
plot(xD_list, dzc(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Center (T%d)', tCenter));
% yline(0,'--','HandleVisibility','off')
xlabel('$x/D$','Interpreter','latex')
ylabel('$\Delta z_c = z_c - z_0$','Interpreter','latex')
legend('Interpreter','latex','Location','best'); box off

yl = ylim;
patch([mergeL mergeR mergeR mergeL],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 0 0], 'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
uistack(findobj(gca,'Type','line'),'top')

%% Optional: ratio metric Rr vs x/D (>=1)
figure('color','white'); hold on
plot(xD_list, Rr(:,1), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Edge (T%d)', tEdge));
plot(xD_list, Rr(:,2), '-o', 'LineWidth', 2, 'DisplayName', sprintf('Center (T%d)', tCenter));
% yline(1,'--','HandleVisibility','off')
xlabel('$x/D$','Interpreter','latex')
ylabel('$R = \max(A_L,A_R)/\min(A_L,A_R)$','Interpreter','latex')
legend('Interpreter','latex','Location','best'); box off

yl = ylim;
patch([mergeL mergeR mergeR mergeL],[yl(1) yl(1) yl(2) yl(2)], ...
      [0 0 0], 'FaceAlpha', 0.06, 'EdgeColor','none', 'HandleVisibility','off');
uistack(findobj(gca,'Type','line'),'top')







%% Functions

function S = wake_asymmetry_metrics(z, d, z0_mode, bounds_mode, varargin)
%WAKE_ASYMMETRY_METRICS  Model-free wake asymmetry metrics from a 1D profile.
%
% Inputs
%   z         : spanwise coordinate (e.g., z/D), column or row
%   d         : deficit signal (e.g., ΔU/U_inf or ΔU/ΔU_max), same size as z
%   z0_mode   : 'peak' (recommended) or 'zero' (use z0=0)
%   bounds_mode: struct with fields:
%       .type  = 'window'
%       .zmin  = lower bound (e.g., -1.5)
%       .zmax  = upper bound (e.g., +1.5)
%
% Name-value options
%   'Threshold' : keep only points where d >= Threshold (default: [])
%   'MinPts'    : minimum points to compute (default: 10)
%
% Outputs (struct S)
%   S.z0      : chosen center for splitting halves
%   S.AL, S.AR: left/right integrated deficit areas
%   S.I       : signed imbalance (AR-AL)/(AR+AL)
%   S.R       : ratio max(AL,AR)/min(AL,AR) (>=1)
%   S.zc      : deficit centroid
%   S.dzc     : zc - z0
%   S.Atot    : total integrated deficit area
%   S.valid   : true if computed successfully

% Defaults
p = inputParser;
p.addParameter('Threshold', [], @(x) isempty(x) || isscalar(x));
p.addParameter('MinPts', 10, @(x) isscalar(x) && x>=3);
p.parse(varargin{:});
thresh = p.Results.Threshold;
minPts = p.Results.MinPts;

S = struct('z0',nan,'AL',nan,'AR',nan,'I',nan,'R',nan,'zc',nan,'dzc',nan,'Atot',nan,'valid',false);

% Vectorize + clean
z = z(:);
d = d(:);
good = isfinite(z) & isfinite(d);
z = z(good);
d = d(good);

if numel(z) < minPts
    return
end

% Sort by z (trapz needs ordered x)
[z, idx] = sort(z);
d = d(idx);

% Apply bounds
if strcmpi(bounds_mode.type,'window')
    inB = (z >= bounds_mode.zmin) & (z <= bounds_mode.zmax);
else
    error('bounds_mode.type not recognized. Use type="window".');
end
z = z(inB);
d = d(inB);

if numel(z) < minPts
    return
end

% Optional threshold (focus on wake core / reduce tail noise)
if ~isempty(thresh)
    keep = d >= thresh;
    z = z(keep);
    d = d(keep);
    if numel(z) < minPts
        return
    end
end

% Choose split location z0
switch lower(z0_mode)
    case 'peak'
        [~, k] = max(d);
        z0 = z(k);
    case 'zero'
        z0 = 0;
    otherwise
        error('z0_mode must be "peak" or "zero".');
end

% Split halves
L = z <= z0;
R = z >= z0;

% If z0 is near one edge, halves may be empty
if sum(L) < 2 || sum(R) < 2
    return
end

AL = trapz(z(L), d(L));
AR = trapz(z(R), d(R));
Atot = AL + AR;

% Guard against tiny totals
if ~isfinite(Atot) || abs(Atot) < eps
    return
end

% Signed imbalance
I = (AR - AL) / (AR + AL);

% Ratio (>=1)
Rratio = max(AL, AR) / max(min(AL, AR), eps);

% Centroid and shift relative to z0
zc = trapz(z, z.*d) / Atot;
dzc = zc - z0;

% Save
S.z0    = z0;
S.AL    = AL;
S.AR    = AR;
S.I     = I;
S.R     = Rratio;
S.Atot  = Atot;
S.zc    = zc;
S.dzc   = dzc;
S.valid = true;
end

