% Attempt to show the far wake can be scaled
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
half_width_factor = 0.54;


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

    % Crop inner side
    raw_wake(Y > 0) = nan;
    massaged_wake(Y > 0) = nan;

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

    % Plot
    contourf(X, Y, massaged_wake)
    plot(x, massaged_bottom_edge, 'color', 'black', 'linewidth', 2)

    % Save to fit a line
    xs(block, :) = x;
    ys(block, :) = massaged_bottom_edge;
    
end


% Fit using all points
x_fit = reshape(xs.', 1, []);
y_fit = reshape(ys.', 1, []);

ft = fittype('(a + b*x) / (1 + c*x)', 'independent','x', 'coefficients',{'a','b','c'});
mdl = fit(x_fit(:), y_fit(:), ft, 'StartPoint',[y(1), 0, 0.1]);

xq = linspace(min(x_fit), max(x_fit), 500);
yq = feval(mdl, xq);

plot(xq, yq, 'k-', 'LineWidth',2, 'DisplayName','fit', 'color', 'red')


hold off
axis equal
set(gca, 'YDir', 'reverse')
ylim([-5, 0])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND TRANSITION FROM 'UNIFORM' TO GAUSSIAN AT EDGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for block = 2:3
   
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y;
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;

    % Delete extra data
    raw_wake(Y > -6) = nan;
    massaged_wake(Y > -6) = nan;

    % Loop through profiles
    for i = 1:length(x)
        % Get profile + centerline velocity
        profile = double(massaged_wake(:, 10));
        centerline = mean_centerline_velocity(block).massaged.u(10);
    
        % Smooth and find where profile recovers
        smooth_profile = sgolayfilt(profile, 2, 3);
        smooth_profile(smooth_profile > centerline) = nan;
        idx = find(~isnan(smooth_profile), 1, 'first');

        % Save
        shear_boundary(block).y(i) = y(idx);
        shear_boundary(block).idx(i) = idx;
    end
end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FAR WAKE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

position_colors = slanCM('thermal-2', 250);

farm_width = 18;

clc; close all
figure('color', 'white')
sgtitle('$21 \leq x \mathbin{/} D \leq 50$', 'interpreter', 'latex')

c = 1;
hold on
for block = 2:3
   
    % Coordinates
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y;
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;

    skip = 2;
    for i = 1:10:length(x)
        Uc = mean_centerline_velocity(block).massaged.u(i);
        scatter(massaged_wake(1:skip:end,i) / Uc, y(1:skip:end) / farm_width, 5, 'filled', 'MarkerFaceColor', position_colors(c, :))
        S = scatter(massaged_wake(shear_boundary(block).idx(i),i) / Uc, y(shear_boundary(block).idx(i)) / farm_width, ...
                    20, 'filled', 'MarkerFaceColor', 'red');
        c = i + 1;
    end

    % Transition region
    max_trans = max(shear_boundary(block).y, [], 'all', 'omitnan');
    min_trans = min(shear_boundary(block).y, [], 'all', 'omitnan');
    max_max(block) = max_trans / farm_width;
    min_min(block) = min_trans / farm_width;
end

% Transition band
x = [0 10];
y1 = min(min_min(2:3));
y2 = max(max_max(2:3));
patch([x(1) x(2) x(2) x(1)], [y1 y1 y2 y2], 'red', ...
      'FaceAlpha', 0.5, 'EdgeColor', 'none');

P = plot([1,1], [-10, 10], 'linewidth', 1, 'color', 'black');
P.Color(4) = 0.5;
uistack(P, 'bottom')
hold off
set(gca, 'YDir', 'reverse')
xlim([0.7, 1.3])
% ylim([-0.7, 0.25])
ylim([-0.75, 0.5])
yline([-0.5, 0.5], 'color', 'black', 'linewidth', 2)
yticks(-0.75:0.25:0.5)
ylabel('$z \mathbin{/} W$', 'interpreter', 'latex', 'FontSize', 16)
xlabel('$\overline{u} \mathbin{/} u_c$', 'Interpreter', 'latex', 'fontsize', 16)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT EDGE OF FAR WAKE PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

farm_width = 18;

figure('color', 'white')

hold on
for block = 2:3
   
    X = data(block).X;
    x = X(1,:);
    Y = data(block).Y;
    y = Y(:,1);

    % Load data
    raw_wake = data(block).u;
    massaged_wake = cleaned(block).u;


    if block == 2
        color = 'red';
    else
        color = 'black';
    end

    for i = 1:5:length(x)
        Uc = mean_centerline_velocity(block).massaged.u(i);
        profile = double(massaged_wake(:, i));
        profile = sgolayfilt(profile, 3, 15);
        local_halfwidth = abs(feval(mdl, x(i)));

        farm_y = (y - shear_boundary(block).y(i)) / farm_width;
        profile(farm_y > 0) = nan;
        profile = (profile - Uc) / u_inf;

        scatter(profile / max(profile, [], 'all', 'omitnan'), (y - shear_boundary(block).y(i)) / local_halfwidth, ...
                10, 'filled');

        % scatter(profile, (y - shear_boundary(block).y(i)) / local_halfwidth, ...
        %         10, 'filled');
    end
end

% Reference gaussian
eta = linspace(-3, 0, 400);             
gauss_ref = exp(-log(2) * eta.^2);      
plot(-gauss_ref + 1.1, eta, 'k-', 'LineWidth', 3, 'DisplayName', 'Gaussian ref')

% yline(turbine_positions / farm_width)

ylim([-3,0])
xlim([-0.1, 1.15])
xticks(0:0.5:1)
yticks(-3:1:0)
hold off
set(gca, 'YDir', 'reverse')

ylabel('$z \mathbin{/} z_{1/2}$', 'interpreter', 'latex', 'FontSize', 16)
xlabel('$\frac{u - u_{c}}{u_{\infty}}$', 'Interpreter', 'latex', 'fontsize', 16)
% xlabel('$\frac{u - u_{c}}{u_{\infty}} \mathbin{/} \left( \frac{u - u_{c}}{u_{\infty}} \right)_{max}$', 'Interpreter', 'latex', 'fontsize', 16)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAR WAKE CENTERLINE VELOCITY (EXTRAPOLATING FAR WAKE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x_fit2 y_fit2 xs ys c

clear x_fit y_fit xs ys c

c = 1;
figure('color', 'white')
hold on
title('Near Wake Extrapolated Amplitude')
for block = 2:3

    if block == 3
        vis = 'on';
    else
        vis = 'off';
    end

    % Plot individual centerlines
    for t = 1:5
        P = plot(mean_centerline_velocity(block).massaged.x, 1 - (wakes(t, block).massaged.centerline_velocity / u_inf), ...
                 'linewidth', 3, 'color', turbine_colors(t,:), ...
                 'DisplayName', sprintf('Turbine %1.0f', t), ...
                 'HandleVisibility', vis);
        P.Color(4) = 1;
    end

    % Plot mean centerline
    plot(mean_centerline_velocity(block).massaged.x, 1 - (mean_centerline_velocity(block).massaged.u / u_inf), ...
        'linewidth', 3, 'HandleVisibility', vis, 'color', 'black', 'DisplayName', 'Mean')
    
    % Save to fit a line
    xs(c, :) = mean_centerline_velocity(block).massaged.x;
    ys(c, :) = 1 - (mean_centerline_velocity(block).massaged.u / u_inf);
    c = c + 1;
end



% Fit using all points
x_fit2 = reshape(xs.', 1, []);
y_fit2 = reshape(ys.', 1, []);


% Jensen normalized deficit shape
jensen_norm = @(k, x) (1 + 2*k.*x).^(-2);
obj = @(k) sum((y_fit2 - jensen_norm(k, x_fit2)).^2);
k0 = 0.07;

% fminsearch is unconstrained; enforce k>0 via reparameterization
% k = exp(p) ensures positivity
obj_p = @(p) sum( (y_fit2 - jensen_norm(exp(p), x_fit2)).^2 );
p0 = log(k0);
p_hat = fminsearch(obj_p, p0);
k_hat = exp(p_hat);

% Best-fit curve on full x range for plotting
yHat = jensen_norm(k_hat, x_fit2);

% Plot
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
plot(x_fit2, yHat, 'color', 'red', 'linestyle', '--',  'LineWidth', 2, 'DisplayName', sprintf('Jensen Optimized\n$k = %1.3f$', k_hat));
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% BPA14 fit (DeltaU/Uinf)_max: optimize k*, epsilon, Ct WITHOUT clipping

xD = x_fit2;        % x/D
y  = y_fit2;   % (DeltaU/Uinf)_max  (dimensionless)

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
xT = x_fit2;        % x̃ = x/D
y  = y_fit2;   % (ΔU/U∞)_max

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
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside');
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\frac{u_{\infty} - u}{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', 16)



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAR WAKE CENTERLINE VELOCITY (FITTING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x_fit y_fit xs ys c

c = 1;
figure('color', 'white')
hold on
title('Far Wake Optimzied Amplitudes')
for block = 2:3

    if block == 3
        vis = 'on';
    else
        vis = 'off';
    end

    % Plot individual centerlines
    for t = 1:5
        P = plot(mean_centerline_velocity(block).massaged.x, 1 - (wakes(t, block).massaged.centerline_velocity / u_inf), ...
                 'linewidth', 3, 'color', turbine_colors(t,:), ...
                 'DisplayName', sprintf('Turbine %1.0f', t), ...
                 'HandleVisibility', vis);
        P.Color(4) = 1;
    end

    % Plot mean centerline
    plot(mean_centerline_velocity(block).massaged.x, 1 - (mean_centerline_velocity(block).massaged.u / u_inf), ...
        'linewidth', 3, 'HandleVisibility', vis, 'color', 'black', 'DisplayName', 'Mean')
    
    % Save to fit a line
    xs(c, :) = mean_centerline_velocity(block).massaged.x;
    ys(c, :) = 1 - (mean_centerline_velocity(block).massaged.u / u_inf);
    c = c + 1;
end



% Fit using all points
x_fit = reshape(xs.', 1, []);
y_fit = reshape(ys.', 1, []);


% Jensen normalized deficit shape
jensen_norm = @(k, x) (1 + 2*k.*x).^(-2);
k_hat = 0.096;
yHat = jensen_norm(k_hat, x_fit);
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
plot(x_fit, yHat, 'color', 'red', 'linestyle', '--',  'LineWidth', 2, 'DisplayName', sprintf('NW Jensen Optimized\n$k = %1.3f$', k_hat));
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% BPA14 fit (DeltaU/Uinf)_max: optimize k*, epsilon, Ct WITHOUT clipping
C_model = @(kstar, eps, Ct, x) 1 - sqrt( 1 - Ct ./ (8*(kstar.*x + eps).^2) );
k_hat = 0.031;
eps_hat = 0.265;
Ct_hat = 0.639;
yHat_full = C_model(k_hat, eps_hat, Ct_hat, x_fit);
bastankhah_label = sprintf('NW Bastankhah Optimized\n$C_t = %1.3f$\n$k = %1.3f$\n$\\varepsilon = %1.3f$', Ct_hat, k_hat, eps_hat);
% Plot
plot(x_fit, yHat_full, 'color', 'green', 'linestyle', '--', 'LineWidth', 2, 'DisplayName', bastankhah_label);
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% Ishihara & Qian-style fit:  A(x̃) = (a + b x̃ + c (1+x̃)^(-2))^(-2)
xT = x_fit;        % x̃ = x/D
y  = y_fit;   % (ΔU/U∞)_max
A_model = @(p, x) (p(1) + p(2).*x + p(3).*(1+x).^(-2)).^(-2);   % p=[a b c]
p_hat(1) = 0.956;
p_hat(2) = 0.196;
p_hat(3) = 0.343;
yHat_full = A_model(p_hat, xT);
% Plot
ishihara_label = sprintf('NW Ishihara Optimized\n$a = %1.3f$\n$b = %1.3f$\n$c = %1.3f$', p_hat(1), p_hat(2), p_hat(3));
plot(xT, yHat_full, 'color', 'magenta', 'linestyle', '--', 'LineWidth', 2, 'DisplayName', ishihara_label);
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% % Exponential fit (all of far wake)
% % Reasonable start guesses
% ft = fittype('A*exp(-k*x)', 'independent','x', 'coefficients',{'A','k'});
% A0 = max(y_fit);
% k0 = 1 / (max(x_fit) - min(x_fit) + eps);
% model = fit(x_fit(:), y_fit(:), ft, ...
%     'StartPoint', [A0, k0], ...
%     'Lower',      [-Inf, 0]);
% xq = linspace(min(xD), max(xD), 400).';
% yq = feval(model, xq);
% % Plot
% exp_label = sprintf('FW Exp fit\n$A = %1.3f$\n$k = %1.3f$', model.A, model.k);
% plot(xq, yq, 'k--', 'LineWidth', 2, 'DisplayName', exp_label)
% plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% % Exponential fit (just block 2)
% x_fit2 = x_fit(1:195);
% y_fit2 = y_fit(1:195);
% ft = fittype('A*exp(-k*x)', 'independent','x', 'coefficients',{'A','k'});
% A0 = max(y_fit2);
% k0 = 1 / (max(x_fit2) - min(x_fit2) + eps);
% model = fit(x_fit2(:), y_fit2(:), ft, ...
%     'StartPoint', [A0, k0], ...
%     'Lower',      [-Inf, 0]);
% xq = linspace(min(xD), max(xD), 400).';
% yq = feval(model, xq);
% % Plot
% exp_label = sprintf('Block 2 fit\n$A = %1.3f$\n$k = %1.3f$', model.A, model.k);
% plot(xq, yq, 'k:', 'LineWidth', 2, 'DisplayName', exp_label)
% plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% % Exponential fit (just block 3)
% x_fit2 = x_fit(196:end);
% y_fit2 = y_fit(196:end);
% ft = fittype('C + A*exp(-k*x)', 'independent','x', ...
%              'coefficients',{'C','A','k'});
% 
% C0 = median(y_fit2(end-20:end), 'omitnan'); % baseline guess from tail
% A0 = max(y_fit2, [], 'omitnan') - C0;
% A0 = max(A0, 1e-6);
% k0 = 0.05;
% 
% model = fit(x_fit2(:), y_fit2(:), ft, ...
%     'StartPoint', [C0, A0, k0], ...
%     'Lower',      [-0.05, 0, 0], ...
%     'Upper',      [ 0.10, 0.20, 1]);
% 
% xq = linspace(min(xD), max(xD), 400).';
% yq = feval(model, xq);
% 
% % Plot
% exp_label = sprintf('Block 3 Exp fit\n$A = %1.3f$\n$k = %1.3f$\n$C = %1.3f$', model.A, model.k, model.C);
% plot(xq, yq, 'k:', 'LineWidth', 2, 'DisplayName', exp_label)



hold off
ylim([-0.05, 0.2])
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside');
% leg.NumColumns = 2;
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\frac{u_{\infty} - u}{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', 16)




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FAR WAKE CENTERLINE VELOCITY (FITTING, EXPONENTIAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear x_fit y_fit xs ys c

c = 1;
figure('color', 'white')
hold on
title('Far Wake Exponential')
for block = 2:3

    if block == 3
        vis = 'on';
    else
        vis = 'off';
    end

    % Plot individual centerlines
    for t = 1:5
        P = plot(mean_centerline_velocity(block).massaged.x, 1 - (wakes(t, block).massaged.centerline_velocity / u_inf), ...
                 'linewidth', 3, 'color', turbine_colors(t,:), ...
                 'DisplayName', sprintf('Turbine %1.0f', t), ...
                 'HandleVisibility', vis);
        P.Color(4) = 1;
    end

    % Plot mean centerline
    plot(mean_centerline_velocity(block).massaged.x, 1 - (mean_centerline_velocity(block).massaged.u / u_inf), ...
        'linewidth', 3, 'HandleVisibility', vis, 'color', 'black', 'DisplayName', 'Mean')
    
    % Save to fit a line
    xs(c, :) = mean_centerline_velocity(block).massaged.x;
    ys(c, :) = 1 - (mean_centerline_velocity(block).massaged.u / u_inf);
    c = c + 1;
end



% Fit using all points
x_fit = reshape(xs.', 1, []);
y_fit = reshape(ys.', 1, []);





% Exponential fit (all of far wake)
% Reasonable start guesses
ft = fittype('A*exp(-k*x)', 'independent','x', 'coefficients',{'A','k'});
A0 = max(y_fit);
k0 = 1 / (max(x_fit) - min(x_fit) + eps);
model = fit(x_fit(:), y_fit(:), ft, ...
    'StartPoint', [A0, k0], ...
    'Lower',      [-Inf, 0]);
xq = linspace(min(xD), max(xD), 400).';
yq = feval(model, xq);
% Plot
exp_label = sprintf('FW Exp fit\n$A = %1.3f$\n$k = %1.3f$', model.A, model.k);
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')
plot(xq, yq, 'g--', 'LineWidth', 2, 'DisplayName', exp_label)
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% Exponential fit (just block 2)
x_fit2 = x_fit(1:195);
y_fit2 = y_fit(1:195);
ft = fittype('A*exp(-k*x)', 'independent','x', 'coefficients',{'A','k'});
A0 = max(y_fit2);
k0 = 1 / (max(x_fit2) - min(x_fit2) + eps);
model = fit(x_fit2(:), y_fit2(:), ft, ...
    'StartPoint', [A0, k0], ...
    'Lower',      [-Inf, 0]);
xq = linspace(min(xD), max(xD), 400).';
yq = feval(model, xq);
% Plot
exp_label = sprintf('Block 2 fit\n$A = %1.3f$\n$k = %1.3f$', model.A, model.k);
plot(xq, yq, 'magenta:', 'LineWidth', 2, 'DisplayName', exp_label)
plot(nan, nan, 'color', 'white', 'DisplayName', ' ')





% Exponential fit (just block 3)
x_fit2 = x_fit(196:end);
y_fit2 = y_fit(196:end);
ft = fittype('C + A*exp(-k*x)', 'independent','x', ...
             'coefficients',{'C','A','k'});

C0 = median(y_fit2(end-20:end), 'omitnan'); % baseline guess from tail
A0 = max(y_fit2, [], 'omitnan') - C0;
A0 = max(A0, 1e-6);
k0 = 0.05;

model = fit(x_fit2(:), y_fit2(:), ft, ...
    'StartPoint', [C0, A0, k0], ...
    'Lower',      [-0.05, 0, 0], ...
    'Upper',      [ 0.10, 0.20, 1]);

xq = linspace(min(xD), max(xD), 400).';
yq = feval(model, xq);

% Plot
exp_label = sprintf('Block 3 Exp fit\n$A = %1.3f$\n$k = %1.3f$\n$C = %1.3f$', model.A, model.k, model.C);
plot(xq, yq, 'r:', 'LineWidth', 2, 'DisplayName', exp_label)



hold off
ylim([-0.05, 0.2])
leg = legend('interpreter', 'latex', 'box', 'off', 'location', 'eastoutside');
% leg.NumColumns = 2;
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', 16)
ylabel('$\frac{u_{\infty} - u}{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', 16)







%% Recovery length per turbine with linear gap fill

thresholds = [0.02, 0.03, 0.04];   % test as many as you want
data_type = 'massaged';

% You need these per turbine:
% x1{t}, y1{t} for window 1 (e.g., 21-30D)
% x2{t}, y2{t} for window 2 (e.g., 41-50D)

% --- USER: define your window bounds (x/D)
win1 = [21, 30];
win2 = [41, 50];

% Example: turbines 1..5
turbines = 1:5;

% Options for gap fill
opts = struct();
opts.extrapolate = false;     % set true if you want to search beyond max(x2)
opts.xMax = 60;               % only used if extrapolate=true
opts.smoothSpan = 1;          % try 5 or 9 if your curves are noisy

% Storage
nT = numel(turbines);
nThr = numel(thresholds);

xRec = nan(nT, nThr);
status = strings(nT, nThr);

% Also store curves (optional)
curves = cell(nT,1);

for ii = 1:nT
    t = turbines(ii);

    % --- USER: provide x1,y1,x2,y2 for this turbine (see extraction helper below)
    % [x1, y1, x2, y2] = get_centerline_deficit_windows(t, win1, win2, data_type);
    x1 = mean_centerline_velocity(2).massaged.x;
    y1 = 1 - wakes(t, 2).(data_type).centerline_velocity/u_inf;
    x2 = mean_centerline_velocity(3).massaged.x;
    y2 = 1 - wakes(t, 3).(data_type).centerline_velocity/u_inf;

    % Choose a common xq grid if you want consistent resolution
    xq = unique([linspace(win1(1), win1(2), 250), linspace(win2(1), win2(2), 250)]).';
    xq = unique([xq; linspace(win1(2), win2(1), 200).']);  % include dense gap

    opts.xq = xq;

    % xq = unique([linspace(21,30,250), linspace(30,41,200), linspace(41,50,250)]).';
    xq = unique([linspace(win1(1),win1(2),250), ...
             linspace(win1(2),win2(1),200), ...
             linspace(win2(1),win2(2),250)]).';


    out = recovery_length_from_windows(x1(:), y1(:), x2(:), y2(:), thresholds, ...
        'xq', xq, ...
        'smoothSpan', 1, ...
        'extrapolate', false);


    xRec(ii,:) = out.xRec(:).';
    status(ii,:) = out.status(:).';

    curves{ii} = out;  %#ok<SAGROW>


    out = recovery_length_from_windows(x1(:), y1(:), x2(:), y2(:), thresholds, ...
    'xq', xq, 'smoothSpan', 1, 'extrapolate', false);

    % Debug print if any NaNs
    if any(isnan(out.xRec))
        fprintf('Turbine %d: NaN recovery for thresholds: ', t);
        disp(thresholds(isnan(out.xRec)).')
    
        % Check window validity
        fprintf('  window1 valid pts = %d, window2 valid pts = %d\n', ...
            sum(isfinite(x1) & isfinite(y1)), sum(isfinite(x2) & isfinite(y2)));
    
        % Show ranges
        fprintf('  y1 range: [%g, %g], y2 range: [%g, %g]\n', ...
            min(y1,[],'omitnan'), max(y1,[],'omitnan'), ...
            min(y2,[],'omitnan'), max(y2,[],'omitnan'));
    
        % Show whether curve ever goes below each threshold
        for kk = 1:numel(thresholds)
            thr = thresholds(kk);
            if all(out.yq > thr)
                fprintf('  Never goes below thr = %.3g (no_cross)\n', thr);
            end
        end
    end

end

% Make a table
T = table(turbines(:), 'VariableNames', {'Turbine'});
for k = 1:nThr
    T.(sprintf('xRec_thr_%g', thresholds(k))) = xRec(:,k);
    T.(sprintf('status_thr_%g', thresholds(k))) = status(:,k);
end
disp(T)



%% Plot recovery distance vs turbine index for each threshold

figure('color','white');
hold on
title('Wake Recovery Length')
for k = 1:nThr
    h = plot(xRec(:,k), turbine_positions, '-o', 'LineWidth', 2, ...
             'DisplayName', sprintf('$\\frac{u_{\\infty} - u}{u_{\\infty}}$ = %.3g', thresholds(k)));
    h.MarkerFaceColor = h.Color;
    plot(nan, nan, 'DisplayName', ' ', 'Color', 'white')
    uistack(h, 'top')
end
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
legend('Location','eastoutside', 'box', 'off', 'interpreter', 'latex'); box off
yticks(turbine_positions)
set(gca, 'YDir', 'reverse')


% Shade regions where we have data
x1 = [1, 21, 41];
x2 = [10, 30, 50];
y = [-15 15];
for i = 1:3
    patch([x1(i) x2(i) x2(i) x1(i)], [y(1) y(1) y(2) y(2)], 'k', ...
          'FaceAlpha', 0.1, 'EdgeColor', 'none', 'handlevisibility', 'off');
end


% yline([-9, 9], 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off')
P = plot([0, 60], [-9, -9], 'color', 'black', 'linewidth', 2, 'HandleVisibility', 'off');
uistack(P, 'bottom')


% Compute the recovery length from the mean centerline velocity
% turbines   = 1:5;                  % adjust
% thresholds = [0, 0.01, 0.05];       % test what you want
% 
% win1 = [21, 30];
% win2 = [41, 50];
% 
% % Collect each turbine's gap-filled curve on this common xq
% Yq = nan(numel(xq), numel(turbines));
% 
% for ii = 1:numel(turbines)
%     t = turbines(ii);
% 
%     x1 = mean_centerline_velocity(2).massaged.x;
%     y1 = 1 - wakes(t, 2).(data_type).centerline_velocity/u_inf;
%     x2 = mean_centerline_velocity(3).massaged.x;
%     y2 = 1 - wakes(t, 3).(data_type).centerline_velocity/u_inf;
% 
%     % Choose a common xq grid if you want consistent resolution
%     xq = unique([linspace(win1(1), win1(2), 250), linspace(win2(1), win2(2), 250)]).';
%     xq = unique([xq; linspace(win1(2), win2(1), 200).']);  % include dense gap
% 
%     out = recovery_length_from_windows(x1(:), y1(:), x2(:), y2(:), thresholds, ...
%         'xq', xq, ...
%         'smoothSpan', 1, ...       % set 5-11 if needed
%         'extrapolate', false);
% 
%     Yq(:,ii) = out.yq;  % store the full gap-filled curve for this turbine
% end
% 
% % Farm-wide mean deficit curve (centerline mean across turbines)
% yFarmMean = mean(Yq, 2, 'omitnan');
% 
% % Now compute recovery distance from the farm-mean curve.
% % Since it's already gap-filled, we can just find threshold crossings directly.
% 
% farmRec = nan(numel(thresholds),1);
% farmStatus = strings(numel(thresholds),1);
% 
% for k = 1:numel(thresholds)
%     thr = thresholds(k);
% 
%     s = yFarmMean - thr;                 % recovered when <= 0
% 
%     if all(~isfinite(s))
%         farmRec(k) = NaN;
%         farmStatus(k) = "no_data";
%         continue
%     end
% 
%     % If already recovered at earliest xq
%     if s(1) <= 0
%         farmRec(k) = xq(1);
%         farmStatus(k) = "lt_min_range";
%         continue
%     end
% 
%     idx = find(s <= 0, 1, 'first');
%     if isempty(idx)
%         farmRec(k) = NaN;
%         farmStatus(k) = "no_cross";
%         continue
%     end
% 
%     % Linear interpolation between idx-1 and idx
%     xL = xq(idx-1); xR = xq(idx);
%     sL = s(idx-1);  sR = s(idx);
%     if (sR - sL) == 0
%         xc = xR;
%     else
%         xc = xL + (0 - sL) * (xR - xL) / (sR - sL);
%     end
%     farmRec(k) = xc;
% 
%     % xline(xc, 'HandleVisibility', 'off')
% 
%     % Classify where it happened
%     if xc >= win1(1) && xc <= win1(2)
%         farmStatus(k) = "measured";
%     elseif xc > win1(2) && xc < win2(1)
%         farmStatus(k) = "gap";
%     elseif xc >= win2(1) && xc <= win2(2)
%         farmStatus(k) = "measured";
%     else
%         farmStatus(k) = "out_of_range";
%     end
% end




axis equal
xlim([0, 53])
ylim([-12, 6])






%% Plot recovery distance vs turbine index for each threshold
% w/ contour

figure('color','white');
hold on
title('Wake Recovery Length')

for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y, cleaned(block).u, ...
             100, 'linestyle', 'none', ...
             'HandleVisibility', 'off')
    clim([4, u_inf])
    colormap(slanCM('gothic'))
end

for k = 1:nThr
    h = plot(xRec(:,k), turbine_positions, '-o', 'LineWidth', 2, ...
             'DisplayName', sprintf('$\\frac{u_{\\infty} - u}{u_{\\infty}}$ = %.3g', thresholds(k)));
    h.MarkerFaceColor = h.Color;
    plot(nan, nan, 'DisplayName', ' ', 'Color', 'white')
    uistack(h, 'top')
end
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
legend('Location','eastoutside', 'box', 'off', 'interpreter', 'latex'); box off
yticks(turbine_positions)
set(gca, 'YDir', 'reverse')



axis equal
xlim([0, 53])
ylim([-12, 6])













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

function out = recovery_length_from_windows(x1, y1, x2, y2, thresholds, varargin)
%RECOVERY_LENGTH_FROM_WINDOWS Estimate recovery distance for multiple thresholds,
% with a linear gap fill between two measurement windows.
%
% Usage:
%   out = recovery_length_from_windows(x1,y1,x2,y2,thresholds)
%   out = recovery_length_from_windows(...,'xq',xq,'xGapFill',xGap,'smoothSpan',5,'extrapolate',true,'xMax',60)

% ---- Parse name-value options ----
p = inputParser;
p.addParameter('xGapFill', [], @(v) isempty(v) || isvector(v));
p.addParameter('xq', [], @(v) isempty(v) || isvector(v));
p.addParameter('extrapolate', false, @(v) islogical(v) && isscalar(v));
p.addParameter('xMax', NaN, @(v) isscalar(v));
p.addParameter('smoothSpan', 1, @(v) isscalar(v) && v >= 1);
p.parse(varargin{:});
opts = p.Results;

% ---- Clean / sort ----
x1 = x1(:); y1 = y1(:);
x2 = x2(:); y2 = y2(:);
thresholds = thresholds(:);

m1 = isfinite(x1) & isfinite(y1);
m2 = isfinite(x2) & isfinite(y2);
x1 = x1(m1); y1 = y1(m1);
x2 = x2(m2); y2 = y2(m2);

[x1, i1] = sort(x1); y1 = y1(i1);
[x2, i2] = sort(x2); y2 = y2(i2);

xEnd1 = max(x1);
xBeg2 = min(x2);

% ---- Defaults ----
if isempty(opts.xGapFill)
    opts.xGapFill = linspace(xEnd1, xBeg2, 200).';
else
    opts.xGapFill = opts.xGapFill(:);
end

if isempty(opts.xq)
    xq = unique([x1; opts.xGapFill; x2]);
else
    xq = unique(opts.xq(:));
end

% ---- Build gap-filled yq ----
yq = nan(size(xq));

% Interp within windows
in1 = (xq >= min(x1)) & (xq <= max(x1));
in2 = (xq >= min(x2)) & (xq <= max(x2));
yq(in1) = interp1(x1, y1, xq(in1), 'linear');
yq(in2) = interp1(x2, y2, xq(in2), 'linear');

% Connect endpoints across the gap
yEnd1 = interp1(x1, y1, xEnd1, 'linear', 'extrap');
yBeg2 = interp1(x2, y2, xBeg2, 'linear', 'extrap');

inGap = (xq > xEnd1) & (xq < xBeg2);
yq(inGap) = yEnd1 + (yBeg2 - yEnd1) .* ((xq(inGap) - xEnd1) ./ (xBeg2 - xEnd1));

% Optional extrapolation (right side)
if opts.extrapolate
    if isnan(opts.xMax), opts.xMax = max(x2); end

    right = (xq > max(x2)) & (xq <= opts.xMax);
    if any(right)
        n = min(5, numel(x2));
        pR = polyfit(x2(end-n+1:end), y2(end-n+1:end), 1);
        yq(right) = polyval(pR, xq(right));
    end
end

% Optional smoothing
if opts.smoothSpan > 1
    yq = smoothdata(yq, 'movmean', opts.smoothSpan, 'omitnan');
end

% ---- Recovery calculations ----
nThr = numel(thresholds);
xRec = nan(nThr,1);
status = strings(nThr,1);
xBracket = nan(nThr,2);

isInMeasured1 = @(x) (x >= min(x1)) && (x <= max(x1));
isInMeasured2 = @(x) (x >= min(x2)) && (x <= max(x2));
isInGapFun    = @(x) (x > xEnd1) && (x < xBeg2);

for k = 1:nThr
    thr = thresholds(k);
    s = yq - thr;                 % recovered when s <= 0

    % If already below threshold at the first available x, recovery happened upstream
    if s(1) <= 0
        xRec(k) = xq(1);
        xBracket(k,:) = [xq(1), xq(1)];
        status(k) = "lt_min_range";  % recovered before our first x
        continue
    end
    
    % Otherwise find first crossing
    idx = find(s <= 0, 1, 'first');
    if isempty(idx)
        status(k) = "no_cross";
        continue
    end


    if idx == 1
        xc = xq(1);
        xL = xq(1); xR = xq(1);
    else
        xL = xq(idx-1); xR = xq(idx);
        sL = s(idx-1);  sR = s(idx);
        if (sR - sL) == 0
            xc = xR;
        else
            xc = xL + (0 - sL) * (xR - xL) / (sR - sL);
        end
    end

    xRec(k) = xc;
    xBracket(k,:) = [xL, xR];

    if isInMeasured1(xc) || isInMeasured2(xc)
        status(k) = "measured";
    elseif isInGapFun(xc)
        status(k) = "gap";
    elseif opts.extrapolate && (xc > max(x2))
        status(k) = "extrapolated";
    else
        status(k) = "gap";
    end
end

out = struct();
out.xq = xq;
out.yq = yq;
out.thresholds = thresholds;
out.xRec = xRec;
out.status = status;
out.xBracket = xBracket;
out.meta = struct('xEnd1',xEnd1,'xBeg2',xBeg2,'yEnd1',yEnd1,'yBeg2',yBeg2);
end

