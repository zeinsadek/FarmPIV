clc; close all; clear

% Paths
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/Natalie_functions')

castaing_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_PDFs/CastaingMatfiles';

% Cases to load
position_filenames = {'z_0', 'z_2', 'z_9', 'z_10'};

% Load matfiles
for i = 1:length(position_filenames)
    disp(i)
    tmp = load(fullfile(castaing_folder, sprintf('CastaingPlotting_%s.mat', position_filenames{i})));
    castaing.(position_filenames{i}) = tmp.output;
end

clear castaing_folder i tmp

%% plot tiles of overlayed x positions per z location

tickFontSize = 8;
labelFontSize = 10;

labels.z_0 = 'Wake: $T_{IV}$';
labels.z_2 = 'Aisle: $T_{III - IV}$';
labels.z_9 = 'Wake: $T_{I}$';
labels.z_10 = 'Aisle: $T_{0 - I}$';

num_xlocations = numel(castaing.(position_filenames{1}).cases);
colors = flipud(slanCM('plasma', num_xlocations));

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 16, 4.5]);
tile = tiledlayout(1, length(position_filenames), 'TileSpacing', 'tight', 'Padding', 'tight');

for i = 1:length(position_filenames)
    tmp = position_filenames{i};
    n_cases = numel(castaing.(tmp).cases);

    h(i) = nexttile;
    set(h(i), 'FontSize', tickFontSize)
    set(h(i), 'TickLabelInterpreter', 'latex');
    hold on
    for j = 1:n_cases
        c = castaing.(tmp).cases(j);

        tmp_x = castaing.(tmp).rho;
        tmp_y = c.lambda_sq;

        idx = find(tmp_y <= 0.01, 1, 'first');
    
        plot(tmp_x, tmp_y, '-o', 'color', colors(j, :), ...
             'markerfacecolor', colors(j, :), 'markersize', 2, 'linewidth', 1, ...
             'HandleVisibility', 'off'); 

        % % Positive region
        % plot(tmp_x(1:idx), tmp_y(1:idx), 'color', colors(j, :), ...
        %      'linewidth', 1, ...
        %      'linestyle', '-', ...
        %      'HandleVisibility', 'off'); 
        % 
        % scatter(tmp_x(1:idx), tmp_y(1:idx), 8, 'o', 'filled', ...
        %         'markerfacecolor', colors(j, :), ...
        %         'handlevisibility', 'off')
        % 
        % 
        % % Negative region (faded)
        % P = plot(tmp_x(idx:end), tmp_y(idx:end), 'color', 'k', ...
        %      'linewidth', 1, ...
        %      'linestyle', '-', ...
        %      'HandleVisibility', 'off'); 
        % P.Color(4) = 0.1;
        % 
        % scatter(tmp_x(idx:end), tmp_y(idx:end), 8, 'o', 'filled', ...
        %         'markerfacecolor', 'k', ...
        %         'handlevisibility', 'off', ...
        %         'MarkerFaceAlpha', 0.1)

    end
    hold off
    box on
    xscale('log')
    yline(0, 'color', 'black', 'LineWidth', 1, 'layer', 'top', 'HandleVisibility', 'off')

    title(labels.(tmp), 'interpreter', 'latex', 'FontSize', labelFontSize)

    yticks(-0.1:0.1:0.6)
    if i ~= 1
        yticklabels({})
    end
end

% x location legend
hold on
for i = 1:num_xlocations
    label = sprintf('$x/D = %1.1f$', castaing.z_0.cases(i).x);
    plot(nan, nan, 'color', colors(i, :), 'displayname', label, ...
         'linewidth', 2)
end
hold off

leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize);
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 15;

linkaxes(h, 'xy')
ylim([-0.05, 0.3])

xlabel(tile, '$\ell \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel(tile, '$\lambda^2$', 'interpreter', 'latex', 'fontsize', labelFontSize)


% Save figure
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/intermittency';
filename = 'SingleFarm_Lam2.pdf';
pause(3)
exportgraphics(fig, fullfile(figure_folder, filename), 'resolution', 600)
close all

















%% single plot with slopes against x/D

z_locations = [0, 1.5, 9, 10.5];

colors = flipud(slanCM('parula', length(z_locations)));

figure('color', 'white')
hold on
for i = 1:length(position_filenames)
    tmp = position_filenames{i};

    % Select cases at the requested z/D
    z_plot = z_locations(i);  % z/D location to plot
    z_tol = 1e-8;
    slope_cases = castaing.(tmp).cases(abs([castaing.(tmp).cases.z] - z_plot) < z_tol);


    % Sort by downstream position
    [x_plot, sort_idx] = sort([slope_cases.x]);
    fits = [slope_cases(sort_idx).fit];
    slope_plot = [fits.slope];

   

    % Plot one marker per measured location
    plot(x_plot, slope_plot, 'o', ...
        'MarkerSize', 6, 'MarkerFaceColor', colors(i, :), ...
        'MarkerEdgeColor', colors(i, :), 'displayname', sprintf('$z/D = %1.1f$', z_plot))

end
hold off
legend('box', 'off', 'interpreter', 'latex', 'location', 'eastoutside')

% xscale('log')


%% Combined plot of tiles and slopes

tickFontSize = 8;
labelFontSize = 10;

x_locations = castaing.(position_filenames{1}).cases.x;
num_xlocations = numel(castaing.(position_filenames{1}).cases);
x_colors = flipud(slanCM('parula', num_xlocations));
z_colors = slanCM('winter', 4);

clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 16, 7])
tiledlayout(2, 4, 'padding', 'compact', 'tilespacing', 'tight');

% Plot lam^2 tiles
for i = 1:length(position_filenames)
    tmp = position_filenames{i};
    n_cases = numel(castaing.(tmp).cases);

    h(i) = nexttile;
    set(h(i), 'FontSize', tickFontSize)
    set(h(i), 'TickLabelInterpreter', 'latex');
    hold on
    for j = 1:n_cases
        c = castaing.(tmp).cases(j);
    
        plot(castaing.(tmp).rho, c.lambda_sq, '-o', 'color', x_colors(j, :), ...
             'markerfacecolor', x_colors(j, :), 'markersize', 2, 'linewidth', 1); 
    end
    hold off
    xscale('log')
end
linkaxes(h(1:4), 'xy')


% Plot slopes vs x
h(5) = nexttile([1, 4]);
set(h(5), 'FontSize', tickFontSize)
set(h(5), 'TickLabelInterpreter', 'latex');
hold on
for i = 1:length(position_filenames)
    tmp = position_filenames{i};

    % Select cases at the requested z/D
    z_plot = z_locations(i);  % z/D location to plot
    z_tol = 1e-8;
    slope_cases = castaing.(tmp).cases(abs([castaing.(tmp).cases.z] - z_plot) < z_tol);


    % Sort by downstream position
    [x_plot, sort_idx] = sort([slope_cases.x]);
    fits = [slope_cases(sort_idx).fit];
    slope_plot = [fits.slope];

  
    % Plot one marker per measured location
    plot(x_plot, slope_plot, 'o', ...
        'MarkerSize', 4, 'MarkerFaceColor', z_colors(i, :), ...
        'linewidth', 2, ...
        'MarkerEdgeColor', z_colors(i, :), 'displayname', sprintf('$z/D = %1.1f$', z_plot))

end
hold off
