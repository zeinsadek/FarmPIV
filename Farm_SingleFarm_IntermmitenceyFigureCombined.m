%%% Plotting intermmittencey figure for single farm paper

clc; close all; clear

% Paths
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/Natalie_functions')

velocity_signal_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/SingleSignalMatfiles';
PDFs_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_PDFs/SinglePDFMatfiles';

% Cases to load
position_filenames = {'x_1_z_0', 'x_1_z_10', 'x_48_z_0'};
position_names = {'Center', 'Edge', 'Far'};

% Load velocity signals
for i = 1:length(position_filenames)
    disp(i)
    tmp = load(fullfile(velocity_signal_folder, sprintf('SingleFarm_VelocitySignal_%s.mat', position_filenames{i})));
    VelocitySignals.(position_names{i}) = tmp.output;
end

% Load velocity signals
for i = 1:length(position_filenames)
    disp(i)
    tmp = load(fullfile(PDFs_folder, sprintf('SingleFarm_PDFs_%s.mat', position_filenames{i})));
    IncrementPDFs.(position_names{i}) = tmp.output;
end

clc;
clear i tmp velocity_signal_folder PDFs_folder


%% Try plotting velocity signals only

clc; close all
figure('color', 'white')

hold on
for i = 1:length(position_filenames)
    tmp_position = position_names{i};
    plot(VelocitySignals.(tmp_position).time_normalized, VelocitySignals.(tmp_position).signal_normalized, ...
         'Displayname', tmp_position, 'linewidth', 2)
end
hold off
xlim([0, 30])
legend('box', 'off', 'location', 'northoutside', 'orientation', 'horizontal')

clear i tmp_position

%% Try plotting PDFs only

tickFontSize = 8;
labelFontSize = 10;

num_rho = length(IncrementPDFs.(position_names{1}));
PDF_colors = flipud(slanCM('parula', num_rho));

clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 7.5, 10])
tiledlayout(length(position_filenames), 1, 'TileSpacing', 'compact', 'Padding', 'tight')

for i = 1:length(position_filenames)
    tmp_position = position_names{i};

    h(i) = nexttile;
    set(h(i), 'FontSize', tickFontSize)
    set(h(i), 'TickLabelInterpreter', 'latex');
    hold on
    for j = 1:length(IncrementPDFs.(tmp_position))
        plot(IncrementPDFs.(tmp_position)(j).xpdf, IncrementPDFs.(tmp_position)(j).pdf, 'linewidth', 2, ...
             'HandleVisibility', 'off', 'color', PDF_colors(j, :))
    end
    hold off
    title(position_names{i}, 'fontsize', labelFontSize, 'interpreter', 'latex')

    yticks(0:0.5:2)
    if i ~= 1
        yticklabels({})
    end
end

% Increment legend
hold on
for i = 1:length(IncrementPDFs.(position_names{1}))
    plot(nan, nan, 'linewidth', 2, ...
         'color', PDF_colors(i, :), ...
         'displayname', num2str(IncrementPDFs.(tmp_position)(i).rho))
end
hold off

linkaxes(h, 'xy')
xlim([-3.5, 3.5])
leg = legend('box', 'off', 'interpreter', 'latex');
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 10;

clear i j h tmp_position leg

%% Combined tiled plot (velocity signal + PDFs)

tickFontSize = 8;
labelFontSize = 10;
linewidth = 1.5;

% Velocity signal colors
signal_colors.Center = '#ffbe0b';
signal_colors.Edge = '#fb5607';
signal_colors.Far = '#8338ec';

% PDF colors
num_rho = length(IncrementPDFs.(position_names{1}));
PDF_colors = flipud(slanCM('plasma', num_rho));


% Coordinate labels
coordinates.Center = '(1, 0)';
coordinates.Edge = '(1, -10.5)';
coordinates.Far = '(48.5, 0)';


clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 15, 11])
tiledlayout(2, length(position_filenames), 'TileSpacing', 'tight', 'Padding', 'compact')

% Velocity signal plot
h(1) = nexttile([1,3]);
hold on
for i = 1:length(position_filenames)
    tmp_position = position_names{i};
    plot(VelocitySignals.(tmp_position).time_normalized, VelocitySignals.(tmp_position).signal_normalized, ...
         'Displayname', tmp_position, ...
         'linewidth', linewidth, ...
         'color', signal_colors.(tmp_position))
end
hold off
box on
grid on

% yline(1, 'color', 'black', 'linewidth', 2, 'layer', 'bottom', 'alpha', 0.25, 'HandleVisibility', 'off')

set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');
xlim([0, 30])
leg = legend('box', 'off', 'location', 'northoutside', 'orientation', 'vertical', 'interpreter', 'latex', 'fontsize', tickFontSize);
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 10;

xlabel('$t \cdot \overline{u} \mathbin{/} D$', 'Interpreter','latex', 'fontsize', labelFontSize)
ylabel('$u \mathbin{/} u_{\infty}$','Interpreter','latex', 'fontsize', labelFontSize)


% PDFs tiles
for i = 1:length(position_filenames)
    tmp_position = position_names{i};

    h(i + 1) = nexttile;
    set(h(i + 1), 'FontSize', tickFontSize)
    set(h(i + 1), 'TickLabelInterpreter', 'latex');
    hold on
    for j = 1:length(IncrementPDFs.(tmp_position))
        plot(IncrementPDFs.(tmp_position)(j).xpdf, IncrementPDFs.(tmp_position)(j).pdf, ...
             'linewidth', linewidth, ...
             'HandleVisibility', 'off', 'color', PDF_colors(j, :))
    end
    hold off
    % title(position_names{i}, 'fontsize', labelFontSize, 'interpreter', 'latex')
    text(-3.2, 1.85, position_names{i}, 'fontsize', tickFontSize, 'interpreter', 'latex')
    text(-3.2, 1.65, coordinates.(tmp_position), 'fontsize', tickFontSize, 'interpreter', 'latex')

    yticks(0:0.5:2)
    if i ~= 1
        yticklabels({})
    end

    if i == 1
        ylabel('PDF', 'fontsize', labelFontSize, 'interpreter', 'latex')
    end

    if i == 2
       xlabel('$\delta u_{\tau}$', 'fontsize', labelFontSize, 'interpreter', 'latex') 
    end
    box on
end

% Increment legend
hold on
for i = 1:length(IncrementPDFs.(position_names{1}))
    plot(nan, nan, 'linewidth', 2, ...
         'color', PDF_colors(i, :), ...
         'displayname', num2str(IncrementPDFs.(tmp_position)(i).rho))
end
hold off

linkaxes(h(2:end), 'xy')
xlim([-3.5, 3.5])
leg = legend('box', 'off', 'interpreter', 'latex', 'fontsize', tickFontSize);
leg.Layout.Tile = 'east';
leg.IconColumnWidth = 10;


clear h i j leg tmp_position



%% Combined tiled plot (velocity signal + PDFs)

tickFontSize = 8;
labelFontSize = 10;
linewidth = 1.25;

% Velocity signal colors
signal_colors.Center = '#02c39a';
signal_colors.Edge = '#0f4c5c';
signal_colors.Far = '#028090';


% PDF colors
num_rho = length(IncrementPDFs.(position_names{1}));
PDF_colors = flipud(slanCM('plasma', num_rho));

% Coordinate labels
coordinates.Center = '(1, 0)';
coordinates.Edge = '(1, -10.5)';
coordinates.Far = '(48.5, 0)';

% Turbine labels
labels.Center = sprintf('Center: $T_{I}$\n$(1, 0)$');
labels.Edge   = sprintf('Edge: $T_{0-I}$\n$(1, -10.5)$');
labels.Far    = sprintf('Far: $T_{I}$\n$(48.5, 0)$');

turbines.Center = '$T_{I}$';
turbines.Edge = '$T_{0-I}$';
turbines.Far = '$T_{I}$';

% For better tile spacing
nCols = length(position_filenames);
topRows = 4;
pdfRows = 5;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 16, 8]);
tiledlayout(topRows + pdfRows, nCols, ...
    'TileSpacing', 'tight', 'Padding', 'compact')


position_names_plottingOrder = {'Center', 'Far', 'Edge'};

% Velocity signal plot
h(1) = nexttile(1, [topRows, nCols]);

hold on
for i = 1:length(position_filenames)
    tmp_position = position_names_plottingOrder{i};

    % Data plot
    label = labels.(tmp_position);
    plot(VelocitySignals.(tmp_position).time_normalized, VelocitySignals.(tmp_position).signal_normalized, ...
         'Displayname', label, ...
         'HandleVisibility', 'off', ...
         'linewidth', linewidth, ...
         'color', signal_colors.(tmp_position))

    % legend entry
    plot(nan, nan, ...
         'Displayname', label, ...
         'linewidth', 2, ...
         'color', signal_colors.(tmp_position))

    plot(nan, nan, 'color', 'white', 'displayname', ' ')

end
hold off
box on
grid on
grid minor

% yline(1, 'color', 'black', 'linewidth', 2, 'layer', 'bottom', 'alpha', 0.25, 'HandleVisibility', 'off')

set(h(1), 'FontSize', tickFontSize)
set(h(1), 'TickLabelInterpreter', 'latex');
xlim([0, 30])
ylim([0, 1.25])
legSignal = legend(h(1), 'Location', 'northeastoutside', ...
    'Orientation', 'vertical', 'Box', 'off', ...
    'Interpreter', 'latex', 'FontSize', tickFontSize);
legSignal.IconColumnWidth = 15;



xlabel('$t \cdot \overline{u} \mathbin{/} D$', 'Interpreter','latex', 'fontsize', labelFontSize)
ylabel('$u \mathbin{/} u_{\infty}$','Interpreter','latex', 'fontsize', labelFontSize)






% PDFs tiles
for i = 1:length(position_filenames)
    tmp_position = position_names{i};

    h(i + 1) = nexttile(topRows * nCols + i, [pdfRows, 1]);
    set(h(i + 1), 'FontSize', tickFontSize)
    set(h(i + 1), 'TickLabelInterpreter', 'latex');
    hold on
    for j = 1:length(IncrementPDFs.(tmp_position))
        plot(IncrementPDFs.(tmp_position)(j).xpdf, smoothdata(IncrementPDFs.(tmp_position)(j).pdf, 'gaussian', 5), ...
             'linewidth', linewidth, ...
             'HandleVisibility', 'off', 'color', PDF_colors(j, :))
    end
    hold off
    text(-3.2, 1.85, sprintf('%s: %s', position_names{i}, turbines.(tmp_position)), 'fontsize', tickFontSize, 'interpreter', 'latex')
    text(-3.2, 1.65, coordinates.(tmp_position), 'fontsize', tickFontSize, 'interpreter', 'latex')

    yticks(0:0.5:2)
    if i ~= 1
        yticklabels({})
    end

    if i == 1
        ylabel('PDF', 'fontsize', labelFontSize, 'interpreter', 'latex')
    end

    if i == 2
       xlabel('$\delta u_{\tau}$', 'fontsize', labelFontSize, 'interpreter', 'latex') 
    end
    box on
end

% Link PDF axes
linkaxes(h(2:end), 'xy')
xlim(h(end), [-3.5, 3.5])

% Discrete colorbar for the PDF curves
rho_values = [IncrementPDFs.(position_names{end}).rho];
nRho = numel(rho_values);

colormap(h(end), PDF_colors(1:nRho, :))
clim(h(end), [0.5, nRho + 0.5])

cb = colorbar(h(end), 'eastoutside');
cb.Direction = 'reverse';
cb.Ticks = 1:nRho;
cb.TickLabels = compose('%g', rho_values);
cb.TickLength = 0;
cb.TickLabelInterpreter = 'latex';
cb.FontSize = tickFontSize;
cb.Label.String = '$\ell \mathbin{/} D$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = labelFontSize;

% Add (a) and (b)
addPanelLabelsFixed(fig, h(1:2), {'a', 'b'}, 'FontSize', 10, 'OffsetPts', [-33,-7])

% Save figure
% figure_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/intermittency';
% filename = 'SingleFarm_Intermittency_Combined.pdf';
% pause(3)
% exportgraphics(fig, fullfile(figure_folder, filename), 'resolution', 600)
% close all

clear h i j tmp_positioncb fig label nCols newHeight pdfRows



%% Functions

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


