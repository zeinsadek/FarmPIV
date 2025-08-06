%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps')
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Data path
experiment = 'Farm2Farm_10D_Gap';
block = 1;

if block == 1
    planes = {'Plane_1_COMBINED', 'Plane_2_COMBINED', 'Plane_3_COMBINED'};
elseif block == 2
    planes = {'Plane_4_COMBINED', 'Plane_5_COMBINED', 'Plane_6_COMBINED'};
elseif block == 3
    planes = {'Plane_7_COMBINED', 'Plane_8_COMBINED', 'Plane_9_COMBINED'};
end

% Paths
mat_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/combined', experiment);
figure_path = fullfile('/Users/zeinsadek/Desktop/Experiments/Farm/figures', experiment);
save_path = fullfile('/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm/blocks', experiment);

if ~exist(save_path, 'dir')
    mkdir(save_path)
end

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
Y = (data(1).Y - 940) / D;

x = X(1,:);
y = Y(:,1);

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CROP DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

components = {'u', 'v', 'uu', 'vv', 'uv'};

% Cutoffs
left_cutoff = 1;
right_cutoff = 4;
top_cutoff = -13;
bottom_cutoff = 4.5;

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
% clear left_idx right_idx top_idx bottom_idx
% clear left_cutoff right_cutoff top_cutoff bottom_cutoff
% clear x y X Y

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

combined.X = combinedX;
combined.Y = combinedY;

clear c component resolution cropped croppedX croppedY

% Shift different planes
offset = 20 * (block - 1);
combined.X = combined.X + offset;


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

component = 'u';

ax = figure();
hold on
title(sprintf("%s: %s", experiment, component), 'interpreter', 'none');

% Plot combiend data
contourf(combined.X, combined.Y, imgaussfilt(combined.(component), 3), 500, 'linestyle', 'none')

if ismember(component, {'v', 'uv'})
    colormap coolwarm
else
    colormap parula
end

% clim([0, 2.5])
% clim([2, 8])

% Plot turbines
color = 'black';
lineWidth = 3;
nacelleLength = 0.25;
for i = 1:5
    % Turbine Positions
    center = -9 + 3 * (i - 1) ;

    % Rotor
    line([0, 0],[center - 0.5, center + 0.5], 'color', color, 'linewidth', lineWidth)

    % Nacelle
    line([0, nacelleLength], [center, center], 'color', color, 'linewidth', lineWidth)

    % Center line
    % yline(center', 'linestyle', '--')
end
hold off
axis equal
colorbar

% Seams
xline(1 + offset)
xline(4 + offset)
xline(7 + offset)
xline(10 + offset)

% Limits
xlim([-1, 51])
ylim([-13, 4.5])

% Ticks
yticks(-12:3:12)
xticks([0,1 + offset:3:12 + offset])

% Labels
fontSize = 16;
xlabel("$x / D$", "interpreter", "latex", "FontSize", fontSize)
ylabel("$z / D$", "interpreter", "latex", "FontSize", fontSize)

set(gca, 'YDir','reverse')

% Save
% file_name = strcat('SingleFarm_Combined_', component, '.png');
% exportgraphics(ax, fullfile(figure_path, file_name), 'resolution', 300);
% close all

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
save(fullfile(save_path, strcat("Block_", num2str(block), "_APPENDED.mat")), 'combined');
fprintf('Data Saved!\n')