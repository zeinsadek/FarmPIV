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


%% Test of violin plot profile

u_inf = 8;

block = 3;
u_profile1 = data(block).u(:,10) / u_inf;
u_profile2 = data(block).u(:,100)/ u_inf;
u_profile3 = data(block).u(:,195)/ u_inf;
y = data(1).Y(:,1);

clc; close all
figure()
hold on
plot(u_profile1, y)
plot(u_profile2, y)
plot(u_profile3, y)
hold off
set(gca, 'YDir', 'reverse')

figure()
hold on
violinplot(1, u_profile1)
violinplot(2, u_profile2)
violinplot(3, u_profile3)
hold off
ylim([0.1, 1.2])

%% Violin plot at specific x-loctions


x_locations = [2.5, 5.5, 8.5, 22.5, 25.5, 28.5, 42.5, 45.5, 48.5] - 0.5;
colors = {'#352473', '#4B209E', '#6909C2', ...
          '#830CCB', '#9A21C9', '#B033C4', ...
          '#C548BA', '#CA67B7', '#C883C5'};

faceAlpha = 0.5;
plot_gap = 1;
labelFontSize = 10;
tickFontSize = 8;

sz = 20;

% Indicies for where the wake centers are and where the aisles are
wake_positions = [-9, -6, -3, 0, 3];
aisle_positions = wake_positions + 1.5;
y = data(block).Y(:,1);

for j = 1:length(wake_positions)
    [~, wake_idx] = min(abs(y - wake_positions(j)));
    [~, aisle_idx] = min(abs(y - aisle_positions(j)));

    wake_indicies(j) = wake_idx;
    aisle_indicies(j) = aisle_idx;
end


clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 7.5, 7])
tiledlayout(1,1,'padding', 'tight')

h(1) = nexttile();
set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');

hold on
for i = 1:length(x_locations)

    x_location = x_locations(i);
    fprintf('x/D = %2.1f\n\n', x_location)

    if x_location >= 1 && x_location <= 10
        block = 1;
    elseif x_location >= 21 && x_location <= 30
        block = 2;
    elseif x_location >= 41 && x_location <= 50
        block = 3;
    end
    disp(block)

    % This is the arbitrary plotting position
    violin_x_positions(i) = i + (plot_gap * (block - 1));

    % Grab data
    y = data(block).Y(:,1);
    x = data(block).X(1,:);

    [~, idx] = min(abs(x - x_location));
    tmp_data = cleaned(block).u;
    tmp_data(data(block).Y < -9) = nan;
    tmp_data(data(block).Y > 3) = nan;
    u_profile = tmp_data(:,idx) / u_inf;

    violinplot(i + (plot_gap * (block - 1)), u_profile, ...
               'FaceColor', colors{i}, ...
               'FaceAlpha', faceAlpha)

    % Plot mean of the entire profile
    % scatter(i + (plot_gap * (block - 1)), mean(u_profile, 'all', 'omitnan'), ...
    %         sz, 'filled', 'MarkerFaceColor', colors{i})

    full_mean_tmp(i) = mean(u_profile, 'all', 'omitnan');

    % Wake-only mean
    % scatter(i + (plot_gap * (block - 1)), mean(u_profile(wake_indicies), 'all', 'omitnan'), ...
    %         sz, 'filled', '^', 'MarkerFaceColor', colors{i})

    wake_mean_tmp(i) = mean(u_profile(wake_indicies), 'all', 'omitnan');

    % Aisle-only mean
    % scatter(i + (plot_gap * (block - 1)), mean(u_profile(aisle_indicies), 'all', 'omitnan'), ...
    %         sz, 'filled', 'square', 'MarkerFaceColor', colors{i})

    aisle_mean_tmp(i) = mean(u_profile(aisle_indicies), 'all', 'omitnan');


end

% Try fitting a polyomial to smooth it out
lineAlpha = 0.25;
lineWidth = 1;
polyOrder = 7;
xfit = linspace(min(violin_x_positions), max(violin_x_positions), 100);

% full_pp = polyfit(violin_x_positions, full_mean_tmp, polyOrder);
% P = plot(xfit, polyval(full_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
% P.Color(4) = lineAlpha;
% uistack(P, 'bottom')
% 
% wake_pp = polyfit(violin_x_positions, wake_mean_tmp, polyOrder);
% P = plot(xfit, polyval(wake_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
% P.Color(4) = lineAlpha;
% uistack(P, 'bottom')
% 
% aisle_pp = polyfit(violin_x_positions, aisle_mean_tmp, polyOrder);
% P = plot(xfit, polyval(aisle_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
% P.Color(4) = lineAlpha;
% uistack(P, 'bottom')

hold off
ylim([0.2, 1.05])
box on

yline(1, 'linewidth', 1)

xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\overline{u} \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

xticks(violin_x_positions)
xticklabels(compose('%.0f', x_locations))


%% STD of profiles plot


clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 4])
tiledlayout(1,1,'padding', 'compact')

h(1) = nexttile();
set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');

hold on

% Storage for all data
x_all = [];
std_all = [];

for block = 1:3
    % Grab data
    y = data(block).Y(:,1);
    x = data(block).X(1,:);

    tmp_data = cleaned(block).u / u_inf;
    tmp_data(data(block).Y < -9) = nan;
    tmp_data(data(block).Y > 3) = nan;

    % Compute STD across y at each x-location
    tmp_std = std(tmp_data, 0, 1, 'omitnan');

    % Plot individual block data
    plot(x, tmp_std, 'color', 'black', 'linewidth', 2)

    % Append to combined arrays
    x_all = [x_all, x];
    std_all = [std_all, tmp_std];
end

% Remove NaNs before fitting
valid = ~isnan(x_all) & ~isnan(std_all);
x_fit_data = x_all(valid);
std_fit_data = std_all(valid);

% Column vectors for fitting
x_fit_data = x_fit_data(:);
std_fit_data = std_fit_data(:);

% Exponential decay model: a*exp(-b*x) + c
exp_model = fittype('a*exp(-b*x) + c', ...
    'independent', 'x', ...
    'coefficients', {'a', 'b', 'c'});

% Initial guesses
a0 = max(std_fit_data) - min(std_fit_data);
b0 = 0.05;
c0 = min(std_fit_data);

fit_opts = fitoptions(exp_model);
fit_opts.StartPoint = [a0, b0, c0];
fit_opts.Lower = [0, 0, 0];   % optional but useful for decay
fit_opts.Upper = [Inf, Inf, Inf];

% Perform fit
[exp_fit, gof] = fit(x_fit_data, std_fit_data, exp_model, fit_opts);

% Plot fitted curve
x_smooth = linspace(min(x_fit_data), max(x_fit_data), 500);
std_smooth = exp_fit.a .* exp(-exp_fit.b .* x_smooth) + exp_fit.c;

plot(x_smooth, std_smooth, '--r', 'linewidth', 2)

% Labels
xlabel('$x/D$', 'interpreter', 'latex')
ylabel('$\sigma_{\overline{u}} \mathbin{/} U_\infty$', 'interpreter', 'latex')

legend({'Block 1', 'Block 2', 'Block 3', 'Exponential fit'}, ...
    'interpreter', 'latex', ...
    'location', 'best')

box on

% Display fit info
disp(exp_fit)
disp(gof)

% Shade the regions where we have data
ax = gca;
yl = ax.YLim;
xmins = [1, 21, 41];
data_shade_color = [0.91, 0.91 ,0.91];

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
xlim([0, 51])




%% Violin plot + inset STD curve

%% Violin plot at specific x-loctions


x_locations = [2.5, 5.5, 8.5, 22.5, 25.5, 28.5, 42.5, 45.5, 48.5] - 0.5;
colors = {'#352473', '#4B209E', '#6909C2', ...
          '#830CCB', '#9A21C9', '#B033C4', ...
          '#C548BA', '#CA67B7', '#C883C5'};

faceAlpha = 0.5;
plot_gap = 1;
labelFontSize = 10;
tickFontSize = 8;

sz = 20;

% Indicies for where the wake centers are and where the aisles are
wake_positions = [-9, -6, -3, 0, 3];
aisle_positions = wake_positions + 1.5;
y = data(block).Y(:,1);

for j = 1:length(wake_positions)
    [~, wake_idx] = min(abs(y - wake_positions(j)));
    [~, aisle_idx] = min(abs(y - aisle_positions(j)));

    wake_indicies(j) = wake_idx;
    aisle_indicies(j) = aisle_idx;
end


clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 7.5, 7])
tiledlayout(1,1,'padding', 'tight')

h(1) = nexttile();
set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');

hold on
for i = 1:length(x_locations)

    x_location = x_locations(i);
    fprintf('x/D = %2.1f\n\n', x_location)

    if x_location >= 1 && x_location <= 10
        block = 1;
    elseif x_location >= 21 && x_location <= 30
        block = 2;
    elseif x_location >= 41 && x_location <= 50
        block = 3;
    end
    disp(block)

    % This is the arbitrary plotting position
    violin_x_positions(i) = i + (plot_gap * (block - 1));

    % Grab data
    y = data(block).Y(:,1);
    x = data(block).X(1,:);

    [~, idx] = min(abs(x - x_location));
    tmp_data = cleaned(block).u;
    tmp_data(data(block).Y < -9) = nan;
    tmp_data(data(block).Y > 3) = nan;
    u_profile = tmp_data(:,idx) / u_inf;

    violinplot(i + (plot_gap * (block - 1)), u_profile, ...
               'FaceColor', colors{i}, ...
               'FaceAlpha', faceAlpha)

    % Plot mean of the entire profile
    scatter(i + (plot_gap * (block - 1)), mean(u_profile, 'all', 'omitnan'), ...
            sz, 'filled', 'MarkerFaceColor', colors{i})

    full_mean_tmp(i) = mean(u_profile, 'all', 'omitnan');

    % Wake-only mean
    scatter(i + (plot_gap * (block - 1)), mean(u_profile(wake_indicies), 'all', 'omitnan'), ...
            sz, 'filled', '^', 'MarkerFaceColor', colors{i})

    wake_mean_tmp(i) = mean(u_profile(wake_indicies), 'all', 'omitnan');

    % Aisle-only mean
    scatter(i + (plot_gap * (block - 1)), mean(u_profile(aisle_indicies), 'all', 'omitnan'), ...
            sz, 'filled', 'square', 'MarkerFaceColor', colors{i})

    aisle_mean_tmp(i) = mean(u_profile(aisle_indicies), 'all', 'omitnan');


end

% Try fitting a polyomial to smooth it out
lineAlpha = 0.25;
lineWidth = 1;
polyOrder = 7;
xfit = linspace(min(violin_x_positions), max(violin_x_positions), 100);

full_pp = polyfit(violin_x_positions, full_mean_tmp, polyOrder);
P = plot(xfit, polyval(full_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
P.Color(4) = lineAlpha;
uistack(P, 'bottom')

wake_pp = polyfit(violin_x_positions, wake_mean_tmp, polyOrder);
P = plot(xfit, polyval(wake_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
P.Color(4) = lineAlpha;
uistack(P, 'bottom')

aisle_pp = polyfit(violin_x_positions, aisle_mean_tmp, polyOrder);
P = plot(xfit, polyval(aisle_pp, xfit), 'color', 'black', 'linewidth', lineWidth);
P.Color(4) = lineAlpha;
uistack(P, 'bottom')

hold off
ylim([0.2, 1.05])
box on

yline(1, 'linewidth', 1)

xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\overline{u} \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)

xticks(violin_x_positions)
xticklabels(compose('%.0f', x_locations))






%%% Inset plot
% horizontal, vertical, height, width
handaxes2 = axes('position', [0.5 0.3 0.45 0.3]);

set(handaxes2, 'FontSize', tickFontSize)
set(handaxes2, 'TickLabelInterpreter', 'latex');

hold on

% Storage for all data
x_all = [];
std_all = [];

for block = 1:3
    % Grab data
    y = data(block).Y(:,1);
    x = data(block).X(1,:);

    tmp_data = cleaned(block).u / u_inf;
    tmp_data(data(block).Y < -9) = nan;
    tmp_data(data(block).Y > 3) = nan;

    % Compute STD across y at each x-location
    tmp_std = std(tmp_data, 0, 1, 'omitnan');

    % Plot individual block data
    plot(x, tmp_std, 'color', 'black', 'linewidth', 2)

    % Append to combined arrays
    x_all = [x_all, x];
    std_all = [std_all, tmp_std];
end

% Remove NaNs before fitting
valid = ~isnan(x_all) & ~isnan(std_all);
x_fit_data = x_all(valid);
std_fit_data = std_all(valid);

% Column vectors for fitting
x_fit_data = x_fit_data(:);
std_fit_data = std_fit_data(:);

% Exponential decay model: a*exp(-b*x) + c
exp_model = fittype('a*exp(-b*x) + c', ...
    'independent', 'x', ...
    'coefficients', {'a', 'b', 'c'});

% Initial guesses
a0 = max(std_fit_data) - min(std_fit_data);
b0 = 0.05;
c0 = min(std_fit_data);

fit_opts = fitoptions(exp_model);
fit_opts.StartPoint = [a0, b0, c0];
fit_opts.Lower = [0, 0, 0];   % optional but useful for decay
fit_opts.Upper = [Inf, Inf, Inf];

% Perform fit
[exp_fit, gof] = fit(x_fit_data, std_fit_data, exp_model, fit_opts);

% Plot fitted curve
x_smooth = linspace(min(x_fit_data), max(x_fit_data), 500);
std_smooth = exp_fit.a .* exp(-exp_fit.b .* x_smooth) + exp_fit.c;

plot(x_smooth, std_smooth, '--r', 'linewidth', 2)

% Labels
xlabel('$x/D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$\sigma_{\overline{u}} \mathbin{/} U_\infty$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% legend({'Block 1', 'Block 2', 'Block 3', 'Exponential fit'}, ...
%     'interpreter', 'latex', ...
%     'location', 'best')

box on

% Display fit info
disp(exp_fit)
disp(gof)

% Shade the regions where we have data
ax = gca;
yl = ax.YLim;
xmins = [1, 21, 41];
data_shade_color = [0.91, 0.91 ,0.91];

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
xlim([0, 51])
