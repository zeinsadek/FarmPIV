%% Farm profiles for single farm paper cartoon

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


%% Where to take profile

x_location_D = 41;
% x_location_D = 21;

if x_location_D >= 1 && x_location_D <= 10
    block = 1;
elseif x_location_D >= 21 && x_location_D <= 30
    block = 2;
elseif x_location_D >= 41 && x_location_D <= 50
    block = 3;
end

x = cleaned(block).X(1,:);
y = cleaned(block).Y(:,1);
u = cleaned(block).u;

[~, idx] = min(abs(x - x_location_D));

nan_mask = ~isnan(u_slice);

u_slice = u(:, idx);
nan_mask = ~isnan(u_slice);

u_slice = u_slice(nan_mask);
y = y(nan_mask);
u_slice = juliaan_smooth(u_slice, 1);
u_slice_sgolay = sgolayfilt(double(u_slice), 2, 27);
% u_slice_sgolay = sgolayfilt(double(u_slice), 2, 57);

clc; close all
figure('color', 'white')
hold on
plot(u_slice_sgolay, y, 'color', 'black', 'linewidth', 4)
% plot(u_slice + 1, y, 'color', 'black', 'linewidth', 4)
hold off
set(gca, 'YDir', 'reverse')
xlim([0, 8.5])
ylim([-13, 6])

