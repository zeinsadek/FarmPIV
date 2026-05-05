%% Plot post-processed hotwire data from Oldenburg

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all; clear
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
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

clear headers opts mean_columns std_columns csv_path csv_dir


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PAPER READY PLOTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Get colors for each downstream plot
palette = 'gem';
colors = slanCM(palette, length(x_positions));

% Plot properties
y_top_limit = 95;
labelFontSize = 10;
tickFontSize = 8;
linewidth = 1;
markersize = 10;

% Shade
x_start = 0;
x_end = 10;
minY = 0.5;
maxY = 1.5;

% Rotor area transparency
rotorAreaAlpha = 0.1;
turbineAlpha = 0.4;

% Legend entry width
legendIconWidth = 10;


% Loop through downstream locations
clc; close all
ax = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 6]);
tile = tiledlayout(1, 2, 'TileSpacing', 'tight');

% Velocity profile
h(1) = nexttile();
set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');
hold on
for l = 1:length(x_positions)

    % Get x-location
    x_position = x_positions(l);
    label = sprintf("$x / D = %2.0f$", x_position * (mesh_spacing / turbine_diameter));

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


% Draw a turbine
turb_x = 5.5;
tower_offset = 0.1;
% Rotor
P = plot([turb_x turb_x], [1.5, 0.5], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;
% Tower
P = plot([turb_x turb_x] + tower_offset, [1, 0], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;
% Nacelle
P = plot([turb_x, turb_x + 2 * tower_offset], [1, 1], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;

% Shading rotor area
P = patch([x_start x_end x_end x_start], [minY minY maxY maxY], ...
      'black', 'FaceAlpha', rotorAreaAlpha, 'EdgeColor', 'none', 'handlevisibility', 'off');
uistack(P, 'bottom')

% Lines at top and bottom tips
yline(1.5, 'linewidth', 0.5, 'color', 'black', 'alpha', 0.15, 'handlevisibility', 'off');
yline(0.5, 'linewidth', 0.5, 'color', 'black', 'alpha', 0.15, 'handlevisibility', 'off');

% Label for shaded region
% plot(nan, nan, 'color', 'white', 'displayname', ' ')
% patch([nan, nan, nan, nan], [nan nan nan nan], ...
%       'black', 'FaceAlpha', rotorAreaAlpha, 'EdgeColor', 'none', 'displayname', 'Rotor Area');

% Horizontal ticks
for i = 2:2:10
    P = yline(i, 'color', 'black', 'linewidth', 0.5, 'alpha', 0.1, 'HandleVisibility', 'off');
    uistack(P, 'bottom')
end

hold off

% Legend
leg = legend('location', 'northwest', 'interpreter', 'latex', 'fontsize', tickFontSize, 'box', 'off');
leg.IconColumnWidth = legendIconWidth;

% Limits and legend
axis square
ylim([0, 12])
xlim([5.2, 10])
yticks(0:2:10)
box on

% Axes labels
xlabel('$\overline{u}$ [m/s]', 'interpreter', 'latex', 'fontsize', labelFontSize)
ylabel('$y \mathbin{/} D$', 'interpreter', 'latex', 'fontsize', labelFontSize)

clear l x_position label prandtl_measurement mean_velocity





% Turbulence intensity
h(2) = nexttile;
set(h(2), 'FontSize', tickFontSize)
set(h(2), 'TickLabelInterpreter', 'latex');
hold on
set(h(2), 'YTickLabel', [])
for l = 1:length(x_positions)

    % Get x-location
    x_position = x_positions(l);
    label = sprintf("$x/D = %2.0f$", x_position * (mesh_spacing / turbine_diameter));

    % Read data
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

% Draw a turbine
turb_x = 0.05;
tower_offset = 0.0024;
% Rotor
P = plot([turb_x turb_x], [1.5, 0.5], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;
% Tower
P = plot([turb_x turb_x] + tower_offset, [1, 0], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;
% Nacelle
P = plot([turb_x, turb_x + 2 * tower_offset], [1, 1], 'color', 'black', 'linewidth', 1, 'handlevisibility', 'off');
P.Color(4) = turbineAlpha;

% Shading
P = patch([x_start x_end x_end x_start], [minY minY maxY maxY], ...
      'black', 'FaceAlpha', rotorAreaAlpha, 'EdgeColor', 'none', 'handlevisibility', 'off');
uistack(P, 'bottom')

% Lines at top and bottom tips
yline(1.5, 'linewidth', 0.5, 'color', 'black', 'alpha', 0.15, 'handlevisibility', 'off');
yline(0.5, 'linewidth', 0.5, 'color', 'black', 'alpha', 0.15, 'handlevisibility', 'off');

% Horizontal ticks
for i = 2:2:10
    P = yline(i, 'color', 'black', 'linewidth', 0.5, 'alpha', 0.1, 'HandleVisibility', 'off');
    uistack(P, 'bottom')
end

hold off


% Limits and legend
axis square
ylim([0, y_top_limit / turbine_diameter])
xlim([0.04, 0.15])
xticks(0.05:0.03:0.2)
box on

% Labels
% xlabel('$Ti$ [\%]', 'interpreter', 'latex', 'fontsize', labelFontSize)
xlabel("$\sqrt{\overline{u' u'}} \mathbin{/} u_{\infty}$", 'interpreter', 'latex', 'fontsize', labelFontSize)
linkaxes(h, 'y')

% Add a) and b) letters
addPanelLabels(h, {'a', 'b'}, 'FontSize', labelFontSize, 'Offset', [-0.1, 1.1])

% Save figure as png and figs
save_path = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper';
clc; fprintf('Saving figures...\n')
pause(3)
exportgraphics(ax, fullfile(save_path, strcat(caze, '_', palette, '.pdf')), 'Resolution', 600);
% savefig(ax, fullfile(save_path, strcat(caze, '.fig')));
close all; clc; fprintf('Figure saved!\n')

clear l label x_position label prandtl_measurement mean_velocity std_velocity turbulence_intensity h 












%% Functions

function addPanelLabels(ax, labels, varargin)
% addPanelLabels(ax, labels) adds (a),(b),... just OUTSIDE top-left of each axes.
% ax     : array of axes handles (e.g., from tiledlayout / findall)
% labels : cellstr like {'a','b','c'} or string array ["a" "b" "c"]
%
% Optional name-value:
% 'Offset'   : [dx dy] in normalized axes units (default [-0.10 1.02])
% 'FontSize' : default 12
% 'FontName' : default 'Times New Roman'

p = inputParser;
addParameter(p,'Offset',[-0.10 1.1]);
addParameter(p,'FontSize', 10);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});
off = p.Results.Offset;

labels = string(labels);
for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Plain parentheses + italic letter:
    % TeX interpreter: \it turns italic ON, \rm returns to roman.
    s = sprintf('(\\ita\\rm)');              % placeholder
    s = sprintf('(\\it%s\\rm)', labels(k));  % actual label

    text(ax(k), off(1), off(2), s, ...
        'Units','normalized', ...
        'Interpreter','tex', ...           % keeps italics control simple
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'Clipping','off');                 % critical: allow outside axes
    end
end


function hAnn = addPanelLabelsFixed(fig, ax, labels, varargin)
% Places (a),(b),... just outside top-left with FIXED physical offsets.
%
% fig    : figure handle (e.g., totalFigure)
% ax     : array of axes handles
% labels : {'a','b',...} or ["a","b",...]
%
% Name-value:
% 'OffsetPts' : [dx dy] in points (default [-10 6])  (left, up)
% 'FontSize'  : default 12
% 'FontName'  : default 'Times New Roman'

p = inputParser;
addParameter(p,'OffsetPts',[-10 6]);
addParameter(p,'FontSize',12);
addParameter(p,'FontName','Times New Roman');
parse(p,varargin{:});

labels = string(labels);
offPts = p.Results.OffsetPts;

ppi = get(0,'ScreenPixelsPerInch');
offPx = offPts/72 * ppi;   % points -> pixels

hAnn = gobjects(numel(ax),1);

for k = 1:numel(ax)
    if ~isgraphics(ax(k),'axes'), continue; end

    % Axes outer position in figure pixels
    op = getpixelposition(ax(k), true);  % [x y w h] in pixels relative to figure

    % Anchor point: outside top-left of axes
    x = op(1) + offPx(1);
    y = op(2) + op(4) + offPx(2);

    % TeX gives roman parentheses + italic letter
    str = sprintf('(\\it%s\\rm)', labels(k));

    hAnn(k) = annotation(fig,'textbox', ...
        'Units','pixels', ...
        'Position',[x y 40 20], ...   % small box; big enough for "(a)"
        'String',str, ...
        'Interpreter','tex', ...
        'FontName',p.Results.FontName, ...
        'FontSize',p.Results.FontSize, ...
        'LineStyle','none', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','bottom');
    end
end



