%% -------------------- LOAD SAVED PLOT DATA --------------------------

clc; close all; clear

matfile_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_PDFs/CastaingMatfiles';
output_file = 'CastaingPlotting_z_10.mat';
loaded = load(fullfile(matfile_folder, output_file), 'output');
output = loaded.output;

%% -------------------- PLOT LAMBDA-SQUARED ---------------------------

style = output.plot;
n_cases = numel(output.cases);
n_rows = ceil(n_cases / style.n_cols);

fig = figure('Color', 'white');
t = tiledlayout(fig, n_rows, style.n_cols);

for j = 1:n_cases
    c = output.cases(j);
    ax = nexttile(t);
    hold(ax, 'on');
    grid(ax, 'on');

    plot(ax, output.rho, c.lambda_sq, '-o', ...
        'Color', c.color, 'MarkerEdgeColor', c.color, ...
        'MarkerFaceColor', c.color, 'MarkerSize', style.marker_size);

    if ~isempty(c.fit.lambda_sq_line)
        plot(ax, c.fit.rho, c.fit.lambda_sq_line, '-', ...
            'Color', style.fit_color);
        text(ax, style.label_xy(2, 1), style.label_xy(2, 2), ...
            sprintf('m = %.5g', c.fit.coefficients(1)));
    end

    ax.XScale = style.x_scale;
    ax.YLim = style.y_lim;
    text(ax, style.label_xy(1, 1), style.label_xy(1, 2), ...
        sprintf('(x,z) = (%g,%g)', c.x, c.z));
end

xlabel(t, style.x_label, 'Interpreter', 'latex');
ylabel(t, style.y_label, 'Interpreter', 'latex');
fontsize(fig, style.font_size, 'points');
fontname(fig, style.font_name);


%% Plot on top of each other

figure('color', 'white')
hold on
for j = 1:n_cases
    c = output.cases(j);

    plot(output.rho, c.lambda_sq, '-o', ...
        'Color', c.color, 'MarkerEdgeColor', c.color, ...
        'MarkerFaceColor', c.color, 'MarkerSize', style.marker_size, 'linewidth', 2);

    if ~isempty(c.fit.lambda_sq_line)
        % plot(c.fit.rho, c.fit.lambda_sq_line, '-', ...
        %     'Color', style.fit_color);
        % text(ax, style.label_xy(2, 1), style.label_xy(2, 2), ...
        %     sprintf('m = %.5g', c.fit.coefficients(1)));
    end

    % ax.XScale = style.x_scale;
    % ax.YLim = style.y_lim;
    % text(ax, style.label_xy(1, 1), style.label_xy(1, 2), ...
    %     sprintf('(x,z) = (%g,%g)', c.x, c.z));
end
hold off

xscale('log')
% xlabel(t, style.x_label, 'Interpreter', 'latex');
% ylabel(t, style.y_label, 'Interpreter', 'latex');
% fontsize(fig, style.font_size, 'points');
% fontname(fig, style.font_name);


%% -------------------- PLOT FIT SLOPE VS. X/D ------------------------

z_plot = 10.5;  % z/D location to plot

% Select cases at the requested z/D
z_tol = 1e-8;
slope_cases = output.cases(abs([output.cases.z] - z_plot) < z_tol);

if isempty(slope_cases)
    error('No saved cases found at z/D = %g.', z_plot);
end

% Sort by downstream position
[x_plot, sort_idx] = sort([slope_cases.x]);
fits = [slope_cases(sort_idx).fit];
slope_plot = [fits.slope];

% Plot one marker per measured location
figure('Color', 'white')
plot(x_plot, slope_plot, 'o', ...
    'MarkerSize', 6, 'MarkerFaceColor', [0.20, 0.35, 0.60], ...
    'MarkerEdgeColor', [0.20, 0.35, 0.60])
grid on

xlabel('$x/D$', 'Interpreter', 'latex')
ylabel('$d\lambda^2/d\log_{10}(\rho/D)$', 'Interpreter', 'latex')
title(sprintf('$z/D = %g$', z_plot), 'Interpreter', 'latex')


