%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/colormaps');
fprintf("All Paths Imported...\n\n");


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data path
mat_path = '/Users/zeinsadek/Downloads/combined';
planes = {'Plane_1_COMBINED', 'Plane_2_COMBINED', 'Plane_3_COMBINED'};

% Figure path
figure_path = '/Users/zeinsadek/Desktop/Experiments/Farm/figures';

for p = 1:length(planes)
    path = fullfile(mat_path, strcat(planes{p}, '.mat'));
    tmp = load(path);
    tmp = tmp.output;
    data(p) = tmp;
end

clear p tmp path


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET COORDINATES (1ST PLANE): CENTERED AND NORMALIZED
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

D = 80;
X = (data(1).X + 200) / D;
Y = (data(1).Y + 940) / D;

x = X(1,:);
y = Y(:,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROP DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

components = {'u', 'v', 'uu', 'vv', 'uv'};

% Cutoffs
left_cutoff = 1;
right_cutoff = 4;
top_cutoff = 13;
bottom_cutoff = -4.5;

% Find closest values in x and y
[~, right_idx] = min(abs(x - right_cutoff));
[~, left_idx] = min(abs(x - left_cutoff));
[~, top_idx] = min(abs(y - top_cutoff));
[~, bottom_idx] = min(abs(y - bottom_cutoff)); 


for p = 1:length(planes)
    for c = 1:length(components)
        tmp = data(p).(components{c});
        cropped(p).(components{c}) = tmp(top_idx:bottom_idx, left_idx:right_idx);
    end
end

croppedX = X(top_idx:bottom_idx, left_idx:right_idx);
croppedY = Y(top_idx:bottom_idx, left_idx:right_idx);

clear p c tmp data
clear left_idx right_idx top_idx bottom_idx
clear left_cutoff right_cutoff top_cutoff bottom_cutoff
clear x y X Y

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c = 1:length(components)
    component = components{c};
    combined.(component) = horzcat(cropped(1).(component), cropped(2).(component), cropped(3).(component));
end

resolution = mean(diff(croppedX(1,:)));
combinedX = horzcat(croppedX, croppedX + range(croppedX(1,:)) + resolution, croppedX + 2*(range(croppedX(1,:)) + resolution));
combinedY = horzcat(croppedY, croppedY, croppedY);

clear c component resolution cropped croppedX croppedY

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


component = 'uv';

ax = figure();
hold on
title(component);

% Plot combiend data
contourf(combinedX, combinedY, combined.(component), 500, 'linestyle', 'none')
colormap parula
% clim([0, 2.5])
clim([-1, 1])

% Plot turbines
color = 'black';
lineWidth = 3;
nacelleLength = 0.25;
for i = 1:5
    % Turbine Positions
    center = -3 + 3 * (i - 1) ;

    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', lineWidth)

    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', lineWidth)
end
hold off
axis equal
colorbar

% Seams
xline(1)
xline(4)
xline(7)
xline(10)

% Limits
xlim([-1, 10])
ylim([-4.5, 13])

% Ticks
yticks(-3:3:12)
xticks([0,1:3:10])

% Labels
fontSize = 16;
xlabel("$x / D$", "interpreter", "latex", "FontSize", fontSize)
ylabel("$y / D$", "interpreter", "latex", "FontSize", fontSize)

% Save
file_name = strcat('SingleFarm_Combined_', component, '.png');
exportgraphics(ax, fullfile(figure_path, file_name), 'resolution', 300);
close all