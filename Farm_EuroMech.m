% Plots for EuroMech Conference

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
experiments = {'SingleFarm', 'Farm2Farm_10D_Gap', 'Farm2Farm_20D_Gap', 'Farm2Farm_40D_Gap'};
blocks = [1,2,3];

% Figure path
figure_folder = '/Users/zeinsadek/Downloads/EuroMech_Figures';

% Turbine/Farm dimensions
D = 80;
Sy = 3;

% Paths
blocks_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks');

% Store all data in structure
for e = 1:length(experiments)
    experiment = experiments{e};
    for b = 1:length(blocks)
        block = blocks(b);
    
        path = fullfile(blocks_path, experiment, sprintf("Block_%1.0f_APPENDED.mat", block));
        tmp = load(path);
        data.(experiment)(block) = tmp.combined;
    
    end
end

clear e p b tmp path experiment block blocks_path


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE A MASSAGED VERSION OF THE DATA AT THE SEAMS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through components
components = {'u', 'v', 'uu', 'vv', 'uv'};

% Loop through cases
for e = 1:length(experiments)
    experiment = experiments{e};

    % Loop through blocks
    for block = 1:3
        for c = 1:length(components)
            component = components{c};
    
            % Pull in block data
            u = data.(experiment)(block).(component);
            X = data.(experiment)(block).X;
            Y = data.(experiment)(block).Y;
    
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
            cleaned.(experiment)(block).(component) = blended;
            cleaned.(experiment)(block).X = X;
            cleaned.(experiment)(block).Y = Y;
        end
    end
end

clear components block c component u X Y blended left right discontinuity 
clear massage_outer_top massage_outer_bottom e


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CONTOURS OF SINGLE CASE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'SingleFarm';
u_inf = 8;

labelFontSize = 10;
tickFontSize = 8;

% Plot contour
clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
tiledlayout(1,1,'padding', 'tight')
nexttile()

hold on
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', tickFontSize)
for block = 1:3
    contourf(cleaned.(experiment)(block).X, cleaned.(experiment)(block).Y, cleaned.(experiment)(block).u / u_inf, ...
             100, 'linestyle', 'none')
end
hold off
axis equal
set(gca, 'YDir', 'reverse')
yticks(-12:3:3)
xlim([-2, 50])
ylim([-12, 4.5])

xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Colorbar
c = colorbar;
c.Label.FontSize = labelFontSize;
c.Label.String = '$u \mathbin{/} u_{\infty}$';
c.Label.Interpreter = 'latex';

c.TickLabelInterpreter = 'latex';
c.FontSize = tickFontSize;
clim([0.3, 1])


% Plane seams
for block = 1:3
    xline([4,7] + 20 * (block - 1), 'Alpha', 0.1)
end


% Plot turbines
turbineLineWidth = 1;
color = 'black';
nacelleLength = 0.3;
for i = 1:5
    % Turbine Positions
    center = -9 + 3 * (i - 1) ;
    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', turbineLineWidth)
    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', turbineLineWidth)
    clear i
end
hold off


% Save figure
exportgraphics(fig, fullfile(figure_folder, 'SingleFarm_u_Contour.pdf'), 'resolution', 600)
clc; close all

clear block c center color C experiment fig nacelleLength turbineLineWidth 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT VERTICAL PROFILES OF ALL CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lw = 1.0;
x_slice_locations = [1, 5, 21, 45];

% %% Pseudo-log spaced profiles
% %%% ChatGPT where to take slices
% % Define total span
% x_min = 1;
% x_max = 50;
% 
% % Create log-spaced positions
% % You can adjust this for more/fewer slices
% n_slices = 5; 
% x_log = logspace(log10(x_min), log10(x_max), n_slices);
% 
% % Round to nearest integers
% x_log_int = unique(round(x_log));
% 
% % Now keep only those values within your available data blocks
% available_blocks = [1:10, 21:30, 41:50];
% x_slices = intersect(x_log_int, available_blocks);
% 
% clc; disp('Suggested x locations for spanwise profiles:')
% disp(x_slices)
% 
% % Take profiles at specific locations (D)
% % x_locations = [1,4,7,10, 21,24,27,30, 41,44,47,50];
% x_locations = x_slices;


% Colors per case
case_colors.SingleFarm = 'black';
case_colors.Farm2Farm_10D_Gap = '#386641';
case_colors.Farm2Farm_20D_Gap = '#6a994e';
case_colors.Farm2Farm_40D_Gap = '#a7c957';

% Case legend names
case_names.SingleFarm = 'Single $\hspace{20pt}$';
case_names.Farm2Farm_10D_Gap = '$10D \hspace{20pt}$';
case_names.Farm2Farm_20D_Gap = '$20D \hspace{20pt}$';
case_names.Farm2Farm_40D_Gap = '$40D \hspace{20pt}$';

fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,18,10]);
tile = tiledlayout(1, length(x_slice_locations), 'padding', 'compact');

for i = 1:length(x_slice_locations)

    x_slice = x_slice_locations(i);

    % Which block to take data from
    if (x_slice >= 1 && x_slice <= 10)
        block = 1;
    elseif (x_slice >= 21 && x_slice <= 30)
        block = 2;
    elseif (x_slice >= 41 && x_slice <= 50)
        block = 3;
    end

    h(i) = nexttile;
    hold on
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca, 'FontSize', tickFontSize)
    for e = 1:length(experiments)
        experiment = experiments{e};

        X = cleaned.(experiment)(block).X;
        Y = cleaned.(experiment)(block).Y;

        x = X(1,:);
        y = Y(:,1);

        % Find index of x slice
        [~, x_idx] = min(abs(x - x_slice));

        u = cleaned.(experiment)(block).u / u_inf;

        % Plot
        plot(u(:, x_idx), y, 'linewidth', lw, 'HandleVisibility', 'off', 'color', case_colors.(experiment))

    end
    hold off
    title(sprintf('$x \\mathbin{/} D = %i$', x_slice), 'interpreter', 'latex')
    set(h(i), 'YDir', 'reverse')
    ylim([-12, 4.5])
    xlim([0, 1])

    if i == 1
        yticks(-12:3:3)
    else
        yticks([])
        ax = gca;
        ax.YAxis.Visible = 'off';
    end
end

% Legend
hold on
for e = 1:length(experiments)
    experiment = experiments{e};
    label = case_names.(experiment);
    plot(nan, nan, 'color', case_colors.(experiment), 'DisplayName', label, 'linewidth', 2)
end
hold off

leg = legend('box', 'off', 'Orientation', 'horizontal', 'interpreter', 'latex');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 20;

linkaxes(h, 'xy')
xlabel(tile, '$u \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'Fontsize', labelFontSize)
ylabel(tile, '$z \mathbin{/} D$', 'interpreter', 'latex', 'Fontsize', labelFontSize)

% Save figure
exportgraphics(fig, fullfile(figure_folder, 'u_Profiles_AllFarms.pdf'), 'resolution', 600)
clc; close all

clear available_blocks ax block e experiment fig h i label leg n_slices tile 
clear u x x_idx x_locations x_log x_log_int x_max x_min x_slice x_slice_locations x_slices X y Y


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CENTER LINE PROFILES OF ALL CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clc; close all
% fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
% tiledlayout(1,1,'padding', 'tight')
% nexttile()
% 
% hold on
% set(gca, 'TickLabelInterpreter', 'latex')
% set(gca, 'FontSize', tickFontSize)
% 
% % Loop through cases
% for e = 1:length(experiments)
%     experiment = experiments{e};
% 
%     % Loop through blocks
%     x_tmp = [];
%     y_tmp = [];
%     for block = 1:3
%         X = cleaned.(experiment)(block).X;
%         Y = cleaned.(experiment)(block).Y;
%         x = X(1,:);
%         y = Y(:,1);
% 
%         % Find index of x slice
%         [~, y_idx] = min(abs(y - (0)));
%         u = cleaned.(experiment)(block).u / u_inf;
% 
%         plot(x, 1-u(y_idx, :), 'linewidth', lw, 'color', case_colors.(experiment))
% 
%         x_tmp = [x_tmp, x];
%         y_tmp = [y_tmp, u(y_idx, :)];
%     end
% 
%     % P = polyfit(x_tmp, y_tmp, 10);
%     % plot(x_tmp, polyval(P, x_tmp))
% end
% 
% hold off
% xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% ylabel('$1 - \frac{u}{u_{\infty}}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
% xlim([0, 50])
% ylim([0, 1])


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT CENTER LINE PROFILES OF ALL CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chat mixed fitting approach

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
tiledlayout(1,1,'padding', 'tight')
nexttile()

hold on
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', tickFontSize)

% Toggle plotting
plot_raw   = true;
plot_fits  = true;

% Loop through cases
for e = 1:length(experiments)
    experiment = experiments{e};
    this_color = case_colors.(experiment);

    % Loop through blocks separately
    for block = 1:3
        X = cleaned.(experiment)(block).X;
        Y = cleaned.(experiment)(block).Y;
        x = X(1,:);
        y = Y(:,1);

        % Find centerline slice
        [~, y_idx] = min(abs(y - 0));
        u = cleaned.(experiment)(block).u / u_inf;
        deficit = 1 - u(y_idx, :);

        % Sort in x just in case
        [x, sort_idx] = sort(x);
        deficit = deficit(sort_idx);

        % Remove duplicate x values if needed
        [x, unique_idx] = unique(x);
        deficit = deficit(unique_idx);

        % Plot raw data
        if plot_raw
            % plot(x, deficit, 'LineWidth', lw, 'Color', this_color, 'HandleVisibility', 'off')
            skip = 10;
            scatter(x(1:skip:end), deficit(1:skip:end), 5, 'filled', 'MarkerFaceColor', this_color, 'HandleVisibility', 'off')
        end

        % Plot fit for this block
        if plot_fits
            x_fit = linspace(min(x), max(x), 300);

            if block == 1
                % Block 1: Exponential fit
                % y = A*exp(-x/L) + C
                exp_fun = fittype('A*exp(-x/L) + C', ...
                    'independent', 'x', 'coefficients', {'A','L','C'});

                A0 = max(deficit) - min(deficit);
                L0 = 3;
                C0 = min(deficit);

                try
                    exp_fit = fit(x(:), deficit(:), exp_fun, ...
                        'StartPoint', [A0, L0, C0], ...
                        'Lower', [0, 0, 0], ...
                        'Upper', [Inf, Inf, 1]);

                    y_fit = exp_fit.A * exp(-x_fit / exp_fit.L) + exp_fit.C;

                    plot(x_fit, y_fit, '-', ...
                        'LineWidth', lw, ...
                        'Color', this_color, ...
                        'HandleVisibility', 'off')

                    fprintf('\n%s block %d exponential fit:\n', experiment, block)
                    fprintf('  A = %.4f, L = %.4f, C = %.4f\n', ...
                        exp_fit.A, exp_fit.L, exp_fit.C)

                catch
                    warning('Exponential fit failed for %s block %d.', experiment, block)
                end

            else
                % Blocks 2 and 3: Linear fit
                % y = m*x + b
                P = polyfit(x, deficit, 1);
                y_fit = polyval(P, x_fit);

                plot(x_fit, y_fit, '-', ...
                    'LineWidth', lw, ...
                    'Color', this_color, ...
                    'HandleVisibility', 'off')

                fprintf('\n%s block %d linear fit:\n', experiment, block)
                fprintf('  slope = %.6f, intercept = %.6f\n', P(1), P(2))
            end
        end
    end
end
hold off

% Legend
hold on
for e = 1:length(experiments)
    experiment = experiments{e};
    label = case_names.(experiment);
    plot(nan, nan, 'color', case_colors.(experiment), 'DisplayName', label, 'linewidth', 2)
end
hold off

leg = legend('box', 'off', 'Orientation', 'horizontal', 'interpreter', 'latex');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 20;


xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$1 - u \mathbin{/} u_{\infty}$', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlim([0, 50])
ylim([0, 1])
box off

% Save figure
% exportgraphics(fig, fullfile(figure_folder, 'u_CenterProfiles_AllFarms.pdf'), 'resolution', 600)
% clc; close all

clear A0 block C0 deficit e exp_fit exp_fun experiment fig label leg log_fit log_fun L0
clear p_spline plot_exp_fit plot_fits plot_log_fit plot_raw plot_spline pp P skip sort_idx Sy
clear this_color u unique_idx x x00 x_fit x_fit_plot x_tmp x_tmp_unique X y y_exp y_fit y_idx
clear y_log y_spline y_tmp y_tmp_unique Y


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ENTRAINMENT (TESTING)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% experiment = 'SingleFarm';
experiment = 'SingleFarm';

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,13,5]);
tiledlayout(1,1,'padding', 'tight')
nexttile()

hold on
set(gca, 'TickLabelInterpreter', 'latex')
set(gca, 'FontSize', tickFontSize)
for block = 1:3

    X = cleaned.(experiment)(block).X;
    Y = cleaned.(experiment)(block).Y;
    x = X(1,:);
    y = Y(:,1);

    u = cleaned.(experiment)(block).u;
    uv = cleaned.(experiment)(block).uv;

    entrainment = uv .* u;
    [e_dx, e_dy] = gradient(entrainment, mean(diff(x)), mean(diff(y)));

    [u_dx, u_dy] = gradient(u, mean(diff(x)), mean(diff(y)));
    advection = u .* u_dx;

    contourf(X, Y, e_dy, 100, 'linestyle', 'none')

end
hold off
axis equal
colorbar()
set(gca, 'YDir', 'reverse')

yticks(-12:3:3)
xlim([-2, 50])
ylim([-12, 4.5])

xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Colorbar
c = colorbar;
c.Label.FontSize = labelFontSize;
c.Label.String = "$\frac{\partial \overline{u'w'} \cdot \overline{u}}{\partial z}$";
c.Label.Interpreter = 'latex';

c.TickLabelInterpreter = 'latex';
c.FontSize = tickFontSize;


% Save figure
% exportgraphics(fig, fullfile(figure_folder, 'Entrainment_SingleFarm_uw_u_dy.pdf'), 'resolution', 600)
% clc; close all


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT ENTRAINMENT PROFILES FOR ALL CASES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


lw = 1.0;
x_slice_locations = [1, 5, 21, 45];

% %% Pseudo-log spaced profiles
% %%% ChatGPT where to take slices
% % Define total span
% x_min = 1;
% x_max = 50;
% 
% % Create log-spaced positions
% % You can adjust this for more/fewer slices
% n_slices = 5; 
% x_log = logspace(log10(x_min), log10(x_max), n_slices);
% 
% % Round to nearest integers
% x_log_int = unique(round(x_log));
% 
% % Now keep only those values within your available data blocks
% available_blocks = [1:10, 21:30, 41:50];
% x_slices = intersect(x_log_int, available_blocks);
% 
% clc; disp('Suggested x locations for spanwise profiles:')
% disp(x_slices)
% 
% % Take profiles at specific locations (D)
% % x_locations = [1,4,7,10, 21,24,27,30, 41,44,47,50];
% x_locations = x_slices;


% Colors per case
case_colors.SingleFarm = 'black';
case_colors.Farm2Farm_10D_Gap = '#386641';
case_colors.Farm2Farm_20D_Gap = '#6a994e';
case_colors.Farm2Farm_40D_Gap = '#a7c957';

% Case legend names
case_names.SingleFarm = 'Single $\hspace{20pt}$';
case_names.Farm2Farm_10D_Gap = '$10D \hspace{20pt}$';
case_names.Farm2Farm_20D_Gap = '$20D \hspace{20pt}$';
case_names.Farm2Farm_40D_Gap = '$40D \hspace{20pt}$';

fig = figure('color', 'white', 'units', 'centimeters', 'position', [10,10,18,10]);
tile = tiledlayout(1, length(x_slice_locations), 'padding', 'compact');

for i = 1:length(x_slice_locations)

    x_slice = x_slice_locations(i);

    % Which block to take data from
    if (x_slice >= 1 && x_slice <= 10)
        block = 1;
    elseif (x_slice >= 21 && x_slice <= 30)
        block = 2;
    elseif (x_slice >= 41 && x_slice <= 50)
        block = 3;
    end

    h(i) = nexttile;
    hold on
    set(gca, 'TickLabelInterpreter', 'latex')
    set(gca, 'FontSize', tickFontSize)
    for e = 1:length(experiments)
        experiment = experiments{e};

        X = cleaned.(experiment)(block).X;
        Y = cleaned.(experiment)(block).Y;

        x = X(1,:);
        y = Y(:,1);

        % Find index of x slice
        [~, x_idx] = min(abs(x - x_slice));

        u = cleaned.(experiment)(block).u;
        uv = cleaned.(experiment)(block).uv;
        entrainment = uv .* u;
        [e_dx, e_dy] = gradient(entrainment, mean(diff(x)), mean(diff(y)));

        [u_dx, u_dy] = gradient(u, mean(diff(x)), mean(diff(y)));
        advection = u .* u_dx;

        % Plot
        plot(-e_dy(:, x_idx), y, 'linewidth', lw, 'HandleVisibility', 'off', 'color', case_colors.(experiment))

    end
    hold off
    title(sprintf('$x \\mathbin{/} D = %i$', x_slice), 'interpreter', 'latex')
    set(h(i), 'YDir', 'reverse')
    ylim([-12, 4.5])
    xlim([-20, 15])

    if i == 1
        yticks(-12:3:3)
    else
        yticks([])
        ax = gca;
        ax.YAxis.Visible = 'off';
    end
end

% Legend
hold on
for e = 1:length(experiments)
    experiment = experiments{e};
    label = case_names.(experiment);
    plot(nan, nan, 'color', case_colors.(experiment), 'DisplayName', label, 'linewidth', 2)
end
hold off

leg = legend('box', 'off', 'Orientation', 'horizontal', 'interpreter', 'latex');
leg.Layout.Tile = 'north';
leg.IconColumnWidth = 20;

linkaxes(h, 'xy')
xlabel(tile, "$- \frac{\partial \, \overline{u'w'} \cdot \overline{u}}{\partial z}$", 'interpreter', 'latex', 'Fontsize', labelFontSize)
ylabel(tile, '$z \mathbin{/} D$', 'interpreter', 'latex', 'Fontsize', labelFontSize)

% Save figure
% exportgraphics(fig, fullfile(figure_folder, 'uv_u_dy_Profiles_AllFarms.pdf'), 'resolution', 600)
% clc; close all

clear available_blocks ax block e experiment fig h i label leg n_slices tile 
clear u x x_idx x_locations x_log x_log_int x_max x_min x_slice x_slice_locations x_slices X y Y

