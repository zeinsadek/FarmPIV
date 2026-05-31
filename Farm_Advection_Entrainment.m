% Compute advection and entrainment terms
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


%% Compute advection terms

levels = 50;

% u * du/dx
figure('color', 'white')
hold on
for block = 1:3
    u = cleaned(block).u;
    X = cleaned(block).X;
    Y = cleaned(block).Y;

    x = X(1,:);
    y = Y(:,1);

    [du_dx, du_dy] = gradient(u, x * 0.08, y * 0.08);

    contourf(X, Y, u .* du_dx, levels, 'linestyle', 'none')
end
hold off
colorbar()
axis equal
colormap(slanCM('coolwarm'))
axis equal
set(gca, 'YDir', 'reverse')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
title('$\overline{u} \cdot \frac{d \overline{u}}{dx}$', 'interpreter', 'latex')


% w * du/dz
figure('color', 'white')
hold on
for block = 1:3
    u = cleaned(block).u;
    v = cleaned(block).v;
    X = cleaned(block).X;
    Y = cleaned(block).Y;

    x = X(1,:);
    y = Y(:,1);

    [du_dx, du_dy] = gradient(u, x * 0.08, y * 0.08);

    contourf(X, Y, v .* du_dy, levels, 'linestyle', 'none')
end
hold off
colorbar()
axis equal
colormap(flipud(slanCM('gothic')))
axis equal
set(gca, 'YDir', 'reverse')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
title('$\overline{w} \cdot \frac{d \overline{u}}{dz}$', 'interpreter', 'latex')


%% Compute entrainment terms

levels = 50;

% u * u'w'
figure('color', 'white')
hold on
for block = 1:3
    u = cleaned(block).u;
    uu = cleaned(block).uu;
    uv = cleaned(block).uv;
    X = cleaned(block).X;
    Y = cleaned(block).Y;

    x = X(1,:);
    y = Y(:,1);

    contourf(X, Y, u .* uv, levels, 'linestyle', 'none')
end
hold off
colorbar()
axis equal
colormap(slanCM('coolwarm'))
axis equal
set(gca, 'YDir', 'reverse')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
title("$\overline{u} \cdot \overline{u' v'}$", 'interpreter', 'latex')


% w * du/dz
figure('color', 'white')
hold on
for block = 1:3
    u = cleaned(block).u;
    uu = cleaned(block).uu;
    v = cleaned(block).v;
    X = cleaned(block).X;
    Y = cleaned(block).Y;

    x = X(1,:);
    y = Y(:,1);

    [du_dx, du_dy] = gradient(u, x * 0.08, y * 0.08);

    contourf(X, Y, u .* uu, levels, 'linestyle', 'none')
end
hold off
colorbar()
axis equal
colormap(flipud(slanCM('gothic')))
axis equal
set(gca, 'YDir', 'reverse')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex')
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex')
title("$\overline{u} \cdot \overline{u' u'}$", 'interpreter', 'latex')