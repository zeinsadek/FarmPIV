%% Plot post-processed hotwire data from Oldenburg

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear
csv_dir = '/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/Hotwire/pp';
csv_name = 'pp-mean-bl_12ms.csv';
csv_path = fullfile(csv_dir, csv_name);
[~, caze, ~] = fileparts(csv_name);

% Load data
opts = detectImportOptions(csv_path);
data = readmatrix(csv_path, opts);
headers = opts.VariableNames;

% Hotwires vertical position [cm]
hotwire_positions = [1.4, 4.4, 7.5, 10.2, 13.7, 16.8, 20.0, 22.8, 25.8, 28.9, 32.1, 35.1, 44.6, 60.5, 74.4, 90];

% Find mean columns
mean_columns = contains(headers, '_mean');
std_columns = contains(headers, '_std');

% Downstream measurement locations [M]
x_positions = data(:,1);

% Strake grid spacing [cm]
mesh_spacing = 14;

% Which columnds to read
mean_indicies = find(mean_columns);
std_indicies = find(std_columns);

% Turbine dimesions [cm]
turbine_diameter = 8;
hub_height = 8;

% Get colors for each downstream plot
colors = (cool(length(x_positions)));

% Plot properties
y_top_limit = 95;
fontsize = 18;
linewidth = 2;
markersize = 50;


clear headers opts mean_columns std_columns csv_path csv_dir


%% Make a 3D plot of where all the measurements were taken

figure()
view(-45, 20)
hold on
for l = 1:length(x_positions)
    x_position = x_positions(l) * mesh_spacing;
    scatter3(ones(size(hotwire_positions)) * x_position, zeros(size(hotwire_positions)), hotwire_positions, ...
             30, 'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', colors(l,:))

    clear l
end
hold off
grid on

% Labels
set(gca, 'YDir', 'reverse')
xlabel('$x$ [cm]', 'interpreter', 'latex', 'fontsize', fontsize)
ylabel('$z$ [cm]', 'interpreter', 'latex', 'fontsize', fontsize)
zlabel('$y$ [cm]', 'interpreter', 'latex', 'fontsize', fontsize)

% Limits
axis equal
xlim([0, 1.1 * (max(x_positions) * mesh_spacing)])
ylim([-300, 300])
zlim([0, 100])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NON-NORMALIZED VELOCITY PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through downstream locations
% figure('color', 'white')
% hold on
% for l = 1:length(x_positions)
% 
%     % Get x-location
%     x_position = x_positions(l);
%     label = sprintf("$x/D = %2.1f$", x_position * (mesh_spacing / turbine_diameter));
% 
%     % Read data
%     mean_velocity = data(l, mean_indicies(2:end-1));
% 
%     % Remove bad point from the tunnel speed we used
%     if strcmp(csv_name, 'pp-mean-bl_12ms.csv')
%         mean_velocity(9) = nan;
%         mean_velocity = fillmissing(mean_velocity, 'spline');
%     end
% 
%     plot(mean_velocity, hotwire_positions / turbine_diameter, 'color', colors(l,:), 'linewidth', linewidth, 'displayname', label)
%     scatter(mean_velocity, hotwire_positions / turbine_diameter, markersize, 'filled', 'markerfacecolor', colors(l,:), 'HandleVisibility', 'off')
% 
% end
% hold off
% 
% % Plot where turbine is
% yline(hub_height / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height + turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height - turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% 
% % Limits and legend
% ylim([0, y_top_limit / turbine_diameter])
% grid on
% legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% % Labels
% xlabel('$u [m/s]$', 'interpreter', 'latex', 'fontsize', fontsize)
% ylabel('$y / D$', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% 
% clear l x_position label prandtl_measurement mean_velocity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NORMALIZED VELOCITY PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through downstream locations
% figure('color', 'white')
% hold on
% for l = 1:length(x_positions)
% 
%     % Get x-location
%     x_position = x_positions(l);
%     label = sprintf("$x/D = %2.1f$", x_position * (mesh_spacing / turbine_diameter));
% 
%     % Read data
%     prandtl_measurement = data(l,mean_indicies(1));
%     mean_velocity = data(l, mean_indicies(2:end-1));
%     normalized_mean_velocity = mean_velocity / prandtl_measurement;
% 
%     % Remove bad point from the tunnel speed we used
%     if strcmp(csv_name, 'pp-mean-bl_12ms.csv')
%         normalized_mean_velocity(9) = nan;
%         normalized_mean_velocity = fillmissing(normalized_mean_velocity, 'spline');
%     end
% 
%     plot(normalized_mean_velocity, hotwire_positions / turbine_diameter, 'color', colors(l,:), 'linewidth', linewidth, 'displayname', label)
%     scatter(normalized_mean_velocity, hotwire_positions / turbine_diameter, markersize, 'filled', 'markerfacecolor', colors(l,:), 'HandleVisibility', 'off')
% 
% end
% hold off
% 
% % Plot where turbine is
% yline(hub_height / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height + turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height - turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% 
% % Limits and legend
% ylim([0, y_top_limit / turbine_diameter])
% grid on
% legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% % Labels
% xlabel('$u / u_{\infty}$', 'interpreter', 'latex', 'fontsize', fontsize)
% ylabel('$y / D$', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% 
% clear l x_position label prandtl_measurement mean_velocity normalized_mean_velocity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TURBULENCE INTENSITY PROFILES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through downstream locations
% figure('color', 'white')
% hold on
% for l = 1:length(x_positions)
% 
%     % Get x-location
%     x_position = x_positions(l);
%     label = sprintf("$x/D = %2.1f$", x_position * (mesh_spacing / turbine_diameter));
% 
%     % Read data
%     % prandtl_measurement = data(l,mean_indicies(1));
%     mean_velocity = data(l, mean_indicies(2:end-1));
%     std_velocity = data(l, std_indicies(2:end-1));
%     turbulence_intensity = std_velocity ./ mean_velocity;
% 
%     % Remove bad point from the tunnel speed we used
%     if strcmp(csv_name, 'pp-mean-bl_12ms.csv')
%         turbulence_intensity(9) = nan;
%         turbulence_intensity = fillmissing(turbulence_intensity, 'spline');
%     end
% 
%     plot(turbulence_intensity, hotwire_positions / turbine_diameter, 'color', colors(l,:), 'linewidth', linewidth, 'displayname', label)
%     scatter(turbulence_intensity, hotwire_positions / turbine_diameter, markersize, 'filled', 'markerfacecolor', colors(l,:), 'HandleVisibility', 'off')
% 
% end
% hold off
% 
% % Plot where turbine is
% yline(hub_height / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height + turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% yline((hub_height - turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
% 
% % Limits and legend
% ylim([0, y_top_limit / turbine_diameter])
% grid on
% legend('location', 'northeast', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% % Labels
% xlabel('$Ti [\%]$', 'interpreter', 'latex', 'fontsize', fontsize)
% ylabel('$y [cm]$', 'interpreter', 'latex', 'fontsize', fontsize)
% 
% clear l label x_position label prandtl_measurement mean_velocity std_velocity turbulence_intensity


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAPER READY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through downstream locations
close all
ax = figure('color', 'white', 'position', [400, 400, 800, 400]);
tiledlayout(1,2)

% Velocity profile
h(1) = nexttile();
hold on
for l = 1:length(x_positions)

    % Get x-location
    x_position = x_positions(l);
    label = sprintf("$x/D = %2.1f$", x_position * (mesh_spacing / turbine_diameter));

    % Read data
    mean_velocity = data(l, mean_indicies(2:end-1));

    % Remove bad point from the tunnel speed we used
    if strcmp(csv_name, 'pp-mean-bl_12ms.csv')
        mean_velocity(9) = nan;
        mean_velocity = fillmissing(mean_velocity, 'spline');
    end

    plot(mean_velocity, hotwire_positions / turbine_diameter, 'color', colors(l,:), 'linewidth', linewidth, 'displayname', label)
    scatter(mean_velocity, hotwire_positions / turbine_diameter, markersize, 'filled', 'markerfacecolor', colors(l,:), 'HandleVisibility', 'off')

end
hold off

% Plot where turbine is
yline(hub_height / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off', 'layer', 'bottom')
yline((hub_height + turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
yline((hub_height - turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')

% Limits and legend
ylim([0, y_top_limit / turbine_diameter])
xlim([4, 10])
grid on
% legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', fontsize)

% Labels
xlabel('$\overline{u} [m/s]$', 'interpreter', 'latex', 'fontsize', fontsize)
ylabel('$y / D$', 'interpreter', 'latex', 'fontsize', fontsize)


clear l x_position label prandtl_measurement mean_velocity


% Turbulence intensity
h(2) = nexttile;
hold on
for l = 1:length(x_positions)

    % Get x-location
    x_position = x_positions(l);
    label = sprintf("$x/D = %2.1f$", x_position * (mesh_spacing / turbine_diameter));

    % Read data
    % prandtl_measurement = data(l,mean_indicies(1));
    mean_velocity = data(l, mean_indicies(2:end-1));
    std_velocity = data(l, std_indicies(2:end-1));
    turbulence_intensity = std_velocity ./ mean_velocity;

    % Remove bad point from the tunnel speed we used
    if strcmp(csv_name, 'pp-mean-bl_12ms.csv')
        turbulence_intensity(9) = nan;
        turbulence_intensity = fillmissing(turbulence_intensity, 'spline');
    end
    
    plot(turbulence_intensity, hotwire_positions / turbine_diameter, 'color', colors(l,:), 'linewidth', linewidth, 'displayname', label)
    scatter(turbulence_intensity, hotwire_positions / turbine_diameter, markersize, 'filled', 'markerfacecolor', colors(l,:), 'HandleVisibility', 'off')

end
hold off

% Plot where turbine is
yline(hub_height / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'HandleVisibility', 'off', 'layer', 'bottom')
yline((hub_height + turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')
yline((hub_height - turbine_diameter/2) / turbine_diameter, 'color', 'black', 'linewidth', linewidth, 'linestyle', '--', 'HandleVisibility', 'off', 'layer', 'bottom')

% Limits and legend
ylim([0, y_top_limit / turbine_diameter])
xlim([0, 0.16])
grid on
legend('location', 'northeastoutside', 'interpreter', 'latex', 'fontsize', fontsize, 'box', 'off')

% Labels
xlabel('$Ti [\%]$', 'interpreter', 'latex', 'fontsize', fontsize)

linkaxes(h, 'y')

clear l label x_position label prandtl_measurement mean_velocity std_velocity turbulence_intensity h 

% Save figure as png and figs
% save_path = '/Users/zeinsadek/Desktop/Experiments/Farm/figures/Inflow';
% clc; fprintf('Saving figures...\n')
% exportgraphics(ax, fullfile(save_path, strcat(caze, '.png')), 'Resolution', 300);
% savefig(ax, fullfile(save_path, strcat(caze, '.fig')));
% close all; clc; fprintf('Figure saved!\n')






