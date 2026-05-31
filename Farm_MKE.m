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


%% Compute MKE terms

levels = 50;

clc;
% MKE terms
for block = 1:3

    u = cleaned(block).u;
    v = cleaned(block).v;

    uu = cleaned(block).uu;
    vv = cleaned(block).vv;
    uv = cleaned(block).uv;

    X = cleaned(block).X;
    Y = cleaned(block).Y;

    x = X(1,:);
    y = Y(:,1);

    % Mean advection
    [du_dx, du_dy] = gradient(u, x * 0.08, y * 0.08);
    mean_advection(block).x = u.^2 .* du_dx;
    mean_advection(block).y = u .* v .* du_dy;


    % Turbulent transport
    U_uu = u .* uu;
    U_uv = u .* uv;

    [U_uu_dx, U_uu_dy] = gradient(U_uu, x * 0.08, y * 0.08);
    [U_uv_dx, U_uv_dy] = gradient(U_uv, x * 0.08, y * 0.08);

    turbulent_transport(block).x = -1 * U_uu_dx;
    turbulent_transport(block).y = -1 * U_uv_dy;


    % Production
    production(block).x = uu .* du_dx;
    production(block).y = uv .* du_dy;

end



%% Plot terms

% Mean advection
figure('color', 'white')
tiledlayout(2,1)
sgtitle('Mean Advection')

h(1) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  mean_advection(block).x, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$u^2 \cdot \partial u / \partial x$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')

h(2) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  mean_advection(block).y, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$u \cdot w \cdot \partial u / \partial z$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')



% Turbulent transport
figure('color', 'white')
tiledlayout(2,1)
sgtitle('Turbulent Transport')

h(1) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  turbulent_transport(block).x, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$-\partial u \cdot \overline{u'u'} / \partial x$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')

h(2) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  turbulent_transport(block).y, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$-\partial u \cdot \overline{u'w'} / \partial z$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')



% Production
figure('color', 'white')
tiledlayout(2,1)
sgtitle('Production')

h(1) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  production(block).x, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$\overline{u'u'} \partial u / \partial x$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')

h(2) = nexttile;
hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y,  production(block).y, levels, 'linestyle', 'none')
end
hold off
axis equal
colorbar()
title("$\overline{u'w'} \partial u / \partial z$", 'interpreter', 'latex')
set(gca, 'YDir', 'reverse')




%% See what profiles look like in the far wake

% Where
block = 1;

% Random idx
idx = 14;

% Quantity
tmp = production(block).y;
y = cleaned(block).Y(:,1);

clc; close all
figure('color', 'white')
plot(tmp(:, idx), y)
yline(-9)
set(gca, 'YDir', 'reverse')



%% Near wake average around each turbine (Test)

block = 1;
turbine_center = -9;
upper_edge = turbine_center - 1.5;
lower_edge = turbine_center + 1.5;

% Normalizing
u_inf = 8;
D = 0.08;
normalization = (u_inf^3) / D;

% Find index of turbine center
X = cleaned(block).X;
Y = cleaned(block).Y;
y = cleaned(block).Y(:,1);

[~, center_idx] = min(abs(y - turbine_center));
[~, upper_idx] = min(abs(y - upper_edge));
[~, lower_edge] = min(abs(y - lower_edge));


% Crop data + coordinates
X = X(upper_idx:lower_edge, :);
Y = Y(upper_idx:lower_edge, :);

cropped_advection_x = mean_advection(block).x(upper_idx:lower_edge, :);
cropped_advection_y = mean_advection(block).y(upper_idx:lower_edge, :);

cropped_transport_x = turbulent_transport(block).x(upper_idx:lower_edge, :);
cropped_transport_y = turbulent_transport(block).y(upper_idx:lower_edge, :);

cropped_produciton_x = production(block).x(upper_idx:lower_edge, :);
cropped_produciton_y = production(block).y(upper_idx:lower_edge, :);


% Average over cropped regions
avg_advection_x = mean(cropped_advection_x, 'all', 'omitnan');
avg_advection_y = mean(cropped_advection_y, 'all', 'omitnan');

avg_transport_x = mean(cropped_transport_x, 'all', 'omitnan');
avg_transport_y = mean(cropped_transport_y, 'all', 'omitnan');

avg_production_x = mean(cropped_produciton_x, 'all', 'omitnan');
avg_production_y = mean(cropped_produciton_y, 'all', 'omitnan');


% Test bar chart
tmp = [avg_advection_x, avg_advection_y; avg_transport_x, avg_transport_y; avg_production_x, avg_production_y];

figure()
bar(tmp, 'stacked')
ylim([-20, 100])






%% Near wake average around each turbine (Looped)

turbine_centers = -9:3:3;
block = 1;

clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 13])
tiledlayout(5,1, 'TileSpacing', 'compact')
sgtitle('Near Wake MKE Budget')

for i = 1:5

    h(i) = nexttile;

    turbine_center = turbine_centers(i);
    upper_edge = turbine_center - 1.5;
    lower_edge = turbine_center + 1.5;
    
    % Normalizing
    u_inf = 8;
    D = 0.08;
    normalization = (u_inf^3) / D;
    
    % Find index of turbine center
    X = cleaned(block).X;
    Y = cleaned(block).Y;
    y = cleaned(block).Y(:,1);
    
    [~, center_idx] = min(abs(y - turbine_center));
    [~, upper_idx] = min(abs(y - upper_edge));
    [~, lower_edge] = min(abs(y - lower_edge));
    
    
    % Crop data + coordinates
    X = X(upper_idx:lower_edge, :);
    Y = Y(upper_idx:lower_edge, :);
    
    cropped_advection_x = mean_advection(block).x(upper_idx:lower_edge, :);
    cropped_advection_y = mean_advection(block).y(upper_idx:lower_edge, :);
    
    cropped_transport_x = turbulent_transport(block).x(upper_idx:lower_edge, :);
    cropped_transport_y = turbulent_transport(block).y(upper_idx:lower_edge, :);
    
    cropped_produciton_x = production(block).x(upper_idx:lower_edge, :);
    cropped_produciton_y = production(block).y(upper_idx:lower_edge, :);
    
    
    % Average over cropped regions
    avg_advection_x = mean(cropped_advection_x, 'all', 'omitnan');
    avg_advection_y = mean(cropped_advection_y, 'all', 'omitnan');
    
    avg_transport_x = mean(cropped_transport_x, 'all', 'omitnan');
    avg_transport_y = mean(cropped_transport_y, 'all', 'omitnan');
    
    avg_production_x = mean(cropped_produciton_x, 'all', 'omitnan');
    avg_production_y = mean(cropped_produciton_y, 'all', 'omitnan');

    % Save to array
    near_wake_MKE(i).advection_x = avg_advection_x;
    near_wake_MKE(i).advection_y = avg_advection_y;

    near_wake_MKE(i).transport_x = avg_transport_x;
    near_wake_MKE(i).transport_y = avg_transport_y;

    near_wake_MKE(i).production_x = avg_production_x;
    near_wake_MKE(i).production_y = avg_production_y;
    
    
    % Test bar chart
    tmp = [avg_advection_x, avg_advection_y; avg_transport_x, avg_transport_y; avg_production_x, avg_production_y];
    

    barh({'$\mathcal{A}$', '$\mathcal{T}$', '$\mathcal{P}$'}, tmp, 'stacked')
    xlim([-12, 65])

    set(gca, 'YTickLabel', {'$\mathcal{A}$', '$\mathcal{T}$', '$\mathcal{P}$'}, 'TickLabelInterpreter', 'latex');
    
    % axis square
    title(sprintf('Turbine %1.0f: $z/D = %1.0f$', i, turbine_center), 'interpreter', 'latex')

    if i < 4
        % yticklabels({})
        xticklabels({})
    end
end

linkaxes(h, 'xy')





%% Farm wake MKE terms


figure('color', 'white')
tiledlayout(1,2)
sgtitle('Far Wake MKE Budget')


turbine_center = -9;
upper_edge = turbine_center - 1.5;
lower_edge = turbine_center + 1.5;

% upper_edge = -9;
% lower_edge = 0;

    
for block = 2:3

    h(i) = nexttile;

    % Normalizing
    u_inf = 8;
    D = 0.08;
    normalization = (u_inf^3) / D;
    
    % Find index of turbine center
    X = cleaned(block).X;
    Y = cleaned(block).Y;
    y = cleaned(block).Y(:,1);
    
    [~, upper_idx] = min(abs(y - upper_edge));
    [~, lower_edge] = min(abs(y - lower_edge));
    
    
    % Crop data + coordinates
    X = X(upper_idx:lower_edge, :);
    Y = Y(upper_idx:lower_edge, :);
    
    cropped_advection_x = mean_advection(block).x(upper_idx:lower_edge, :);
    cropped_advection_y = mean_advection(block).y(upper_idx:lower_edge, :);
    
    cropped_transport_x = turbulent_transport(block).x(upper_idx:lower_edge, :);
    cropped_transport_y = turbulent_transport(block).y(upper_idx:lower_edge, :);
    
    cropped_produciton_x = production(block).x(upper_idx:lower_edge, :);
    cropped_produciton_y = production(block).y(upper_idx:lower_edge, :);
    
    
    % Average over cropped regions
    avg_advection_x = mean(cropped_advection_x, 'all', 'omitnan');
    avg_advection_y = mean(cropped_advection_y, 'all', 'omitnan');
    
    avg_transport_x = mean(cropped_transport_x, 'all', 'omitnan');
    avg_transport_y = mean(cropped_transport_y, 'all', 'omitnan');
    
    avg_production_x = mean(cropped_produciton_x, 'all', 'omitnan');
    avg_production_y = mean(cropped_produciton_y, 'all', 'omitnan');

     % Save to array
    far_wake_MKE(block).advection_x = avg_advection_x;
    far_wake_MKE(block).advection_y = avg_advection_y;

    far_wake_MKE(block).transport_x = avg_transport_x;
    far_wake_MKE(block).transport_y = avg_transport_y;

    far_wake_MKE(block).production_x = avg_production_x;
    far_wake_MKE(block).production_y = avg_production_y;
    
    
    % Test bar chart
    tmp = [avg_advection_x, avg_advection_y; avg_transport_x, avg_transport_y; avg_production_x, avg_production_y];
    

    bar({'Advection', 'Transport', 'Production'}, tmp, 'stacked')
    ylim([-20, 100])

    axis square
    title(sprintf('Block %1.0f', block), 'interpreter', 'latex')

    if block > 2
        yticklabels({})
    end


end





%% Near wake average around each turbine: Plotted by term

turbine_centers = -9:3:3;

clear h

clc; close all
figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 10])
tiledlayout(3,2, 'TileSpacing', 'compact')
% sgtitle('Near Wake MKE Budget')

% turbine_names = {'T1', 'T2', 'T3', 'T4', 'T5'};
turbine_names = {'$T_{\mathrm{I}$', ...
                 '$T_{\mathrm{II}$', ...
                 '$T_{\mathrm{III}$', ...
                 '$T_{\mathrm{IV}$', ...
                 '$T_{\mathrm{V}$'};

% Different colors
advection_x_color = '#005f60';
advection_y_color = '#fd5901';

transport_x_color = '#0a2344';
transport_y_color = '#fd0363';

production_x_color = '#249ea0';
production_y_color = '#faab36';


% Make array for data being plotted
near_advection_tmp = nan(5,2);
near_transport_tmp = nan(5,2);
near_production_tmp = nan(5,2);

for i = 1:5
    near_advection_tmp(i,1) = near_wake_MKE(i).advection_x;
    near_advection_tmp(i,2) = near_wake_MKE(i).advection_y;

    near_transport_tmp(i,1) = near_wake_MKE(i).transport_x;
    near_transport_tmp(i,2) = near_wake_MKE(i).transport_y;

    near_production_tmp(i,1) = near_wake_MKE(i).production_x;
    near_production_tmp(i,2) = near_wake_MKE(i).production_y;
end


%%% NEAR WAKE
% Advection
h(1) = nexttile(1);
hold on
bh = barh(turbine_names, near_advection_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(advection_x_color);
bh(2).CData = hex2rgb(advection_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
ylabel('Advection', 'interpreter', 'latex')

% Title
title({'Individual wakes', '$x \mathbin{/}D = 1 \, - \, 10$'}, 'interpreter', 'latex')




% Production
h(3) = nexttile(3);
hold on
bh = barh(turbine_names, near_production_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(production_x_color);
bh(2).CData = hex2rgb(production_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
ylabel('Production', 'interpreter', 'latex')



% Transport
h(2) = nexttile(5);
hold on
bh = barh(turbine_names, near_transport_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')
bh(1).ShowBaseLine = 'off';

% Colors
bh(1).CData = hex2rgb(transport_x_color);
bh(2).CData = hex2rgb(transport_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
xlim([-12, 12])
ylabel('Transport', 'interpreter', 'latex')





%%% FAR WAKE
far_names = {'Block 2', 'Block 3'};

% Make array for data being plotted
far_advection_tmp = nan(2,2);
far_transport_tmp = nan(2,2);
far_production_tmp = nan(2,2);

for i = 2:3
    far_advection_tmp(i - 1,1) = far_wake_MKE(i).advection_x;
    far_advection_tmp(i - 1,2) = far_wake_MKE(i).advection_y;

    far_transport_tmp(i - 1,1) = far_wake_MKE(i).transport_x;
    far_transport_tmp(i - 1,2) = far_wake_MKE(i).transport_y;

    far_production_tmp(i - 1,1) = far_wake_MKE(i).production_x;
    far_production_tmp(i - 1,2) = far_wake_MKE(i).production_y;
end


% Advection
h(4) = nexttile(2);
hold on
bh = barh(far_names, far_advection_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(advection_x_color);
bh(2).CData = hex2rgb(advection_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', far_names, 'TickLabelInterpreter', 'latex');

% Legend
leg = legend({'$\mathcal{A}_{x}$', '$\mathcal{A}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'southeast');
leg.IconColumnWidth = 7;

% Title
title({'Farm wake edge', '$x \mathbin{/}D = 21 \, - \, 50$'}, 'interpreter', 'latex')



% Production
h(5) = nexttile(4);
hold on
bh = barh(far_names, far_production_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(production_x_color);
bh(2).CData = hex2rgb(production_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', far_names, 'TickLabelInterpreter', 'latex');

% Legend
leg = legend({'$\mathcal{P}_{x}$', '$\mathcal{P}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'east');
leg.IconColumnWidth = 7;





% Transport
h(6) = nexttile(6);
hold on
bh = barh(far_names, far_transport_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')
bh(1).ShowBaseLine = 'off';

% Colors
bh(1).CData = hex2rgb(transport_x_color);
bh(2).CData = hex2rgb(transport_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', far_names, 'TickLabelInterpreter', 'latex');

% Legend
leg = legend({'$\mathcal{T}_{x}$', '$\mathcal{T}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'east');
leg.IconColumnWidth = 7;





linkaxes([h(1), h(3)], 'x')
xlim([0, 70])

linkaxes([h(4), h(5)], 'x')
xlim([0, 26])

linkaxes([h(2), h(6)], 'x')
xlim([-12, 12])





%% Near wake average around each turbine: Plotted by term NEW 3x3 FORMAT
%%% THIS ONE

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 13, 9]);
tile = tiledlayout(3,3, 'TileSpacing', 'tight', 'padding', 'tight');

% Turbine names in roman numerals
turbine_names = {'$T_{\mathrm{I}}$', ...
                 '$T_{\mathrm{II}}$', ...
                 '$T_{\mathrm{III}}$', ...
                 '$T_{\mathrm{IV}}$', ...
                 '$T_{\mathrm{V}}$'};

% Different colors
advection_x_color = '#005f60';
advection_y_color = '#fd5901';

transport_x_color = '#0a2344';
transport_y_color = '#fd0363';

production_x_color = '#249ea0';
production_y_color = '#faab36';

% Fontsizes
tickFontSize = 8;
labelFontSize = 10;
legendIconWidth = 9;


% Make array for data being plotted
near_advection_tmp = nan(5,2);
near_transport_tmp = nan(5,2);
near_production_tmp = nan(5,2);

for i = 1:5
    near_advection_tmp(i,1) = near_wake_MKE(i).advection_x;
    near_advection_tmp(i,2) = near_wake_MKE(i).advection_y;

    near_transport_tmp(i,1) = near_wake_MKE(i).transport_x;
    near_transport_tmp(i,2) = near_wake_MKE(i).transport_y;

    near_production_tmp(i,1) = near_wake_MKE(i).production_x;
    near_production_tmp(i,2) = near_wake_MKE(i).production_y;
end


% Make array for data being plotted
far_advection_tmp = nan(2,2);
far_transport_tmp = nan(2,2);
far_production_tmp = nan(2,2);

for i = 2:3
    far_advection_tmp(i - 1,1) = far_wake_MKE(i).advection_x;
    far_advection_tmp(i - 1,2) = far_wake_MKE(i).advection_y;

    far_transport_tmp(i - 1,1) = far_wake_MKE(i).transport_x;
    far_transport_tmp(i - 1,2) = far_wake_MKE(i).transport_y;

    far_production_tmp(i - 1,1) = far_wake_MKE(i).production_x;
    far_production_tmp(i - 1,2) = far_wake_MKE(i).production_y;
end






%%% NEAR WAKE: Advection
h(1) = nexttile(1);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(turbine_names, near_advection_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(advection_x_color);
bh(2).CData = hex2rgb(advection_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
ylabel('Advection xx', 'interpreter', 'latex', 'fontsize', labelFontSize)
box on

% Title
title('$1 \leq x \mathbin{/}D \leq 10$', 'interpreter', 'latex', 'fontsize', labelFontSize)



%%% Block 2: Advection
block = 2;
h(2) = nexttile(2);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).advection_x, far_wake_MKE(block).advection_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(advection_x_color);
bh(2).CData = hex2rgb(advection_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
box on

% Title
title('$21 \leq x \mathbin{/}D \leq 30$', 'interpreter', 'latex', 'fontsize', labelFontSize)



%%% Block 3: Advection
block = 3;
h(3) = nexttile(3);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).advection_x, far_wake_MKE(block).advection_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(advection_x_color);
bh(2).CData = hex2rgb(advection_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
box on

% Title
title('$41 \leq x \mathbin{/}D \leq 50$', 'interpreter', 'latex', 'fontsize', labelFontSize)

% Legend
leg = legend({'$\mathcal{A}_{x}$', '$\mathcal{A}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'east', 'fontsize', labelFontSize);
leg.IconColumnWidth = legendIconWidth;


%%% Link xlims for advection
linkaxes([h(1), h(2), h(3)], 'x')
xlim([0, 70])





%%% NEAR WAKE: Production
h(4) = nexttile(4);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(turbine_names, near_production_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(production_x_color);
bh(2).CData = hex2rgb(production_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
ylabel('Production xx', 'interpreter', 'latex', 'fontsize', labelFontSize)
box on




%%% Block 2: Production
block = 2;
h(5) = nexttile(5);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).production_x, far_wake_MKE(block).production_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(production_x_color);
bh(2).CData = hex2rgb(production_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
box on



%%% Block 3: Production
block = 3;
h(6) = nexttile(6);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).production_x, far_wake_MKE(block).production_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(production_x_color);
bh(2).CData = hex2rgb(production_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
box on

% Legend
leg = legend({'$\mathcal{P}_{x}$', '$\mathcal{P}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'east', 'fontsize', labelFontSize);
leg.IconColumnWidth = legendIconWidth;


%%% Link xlims for advection
linkaxes([h(4), h(5), h(6)], 'x')
xlim([0, 9])






%%% NEAR WAKE: Transport
h(7) = nexttile(7);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(turbine_names, near_transport_tmp, 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(transport_x_color);
bh(2).CData = hex2rgb(transport_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', turbine_names, 'TickLabelInterpreter', 'latex');
ylabel('Transport xx', 'interpreter', 'latex', 'fontsize', labelFontSize)
bh(1).ShowBaseLine = 'off';
box on




%%% Block 2: Transport
block = 2;
h(8) = nexttile(8);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).transport_x, far_wake_MKE(block).transport_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(transport_x_color);
bh(2).CData = hex2rgb(transport_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
bh(1).ShowBaseLine = 'off';
box on

% Xlabel for everything
xlabel(tile, '$[m^2 \mathbin{/} s^3]$', 'interpreter', 'latex', 'fontsize', labelFontSize)



%%% Block 3: Transport
block = 3;
h(9) = nexttile(9);
set(gca, 'TickLabelInterpreter', 'latex');
ax = gca;
ax.FontSize = tickFontSize;

hold on
bh = barh(1, [far_wake_MKE(block).transport_x, far_wake_MKE(block).transport_y], 'stacked');
set(bh, 'EdgeColor', 'none');
set(bh, 'FaceColor', 'Flat')

% Colors
bh(1).CData = hex2rgb(transport_x_color);
bh(2).CData = hex2rgb(transport_y_color);
hold off
set(gca, 'YDir', 'reverse')
set(gca,'ytick',[])
bh(1).ShowBaseLine = 'off';
box on

% Legend
leg = legend({'$\mathcal{T}_{x}$', '$\mathcal{T}_{z}$'}, 'box', 'off', ...
             'orientation', 'vertical', 'interpreter', 'latex', ...
             'location', 'east', 'fontsize', labelFontSize);
leg.IconColumnWidth = legendIconWidth;


%%% Link xlims for advection
linkaxes([h(7), h(8), h(9)], 'x')
xlim([-14, 14])



% Save figure
clc; 
save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MKE';
figure_name = 'SingleFarmPaper_MKE_Terms.pdf';
pause(3)
fprintf('Saving figure...\n')
exportgraphics(fig, fullfile(save_folder, figure_name), 'resolution', 600)
clc; fprintf('Figure Saved!\n')







%% Contour cartoon for where we draw our control volumes

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 12.58, 2.5]);
tile = tiledlayout(1,1,'padding', 'tight');

nexttile()
ax = gca;
ax.FontSize = tickFontSize;
set(gca, 'TickLabelInterpreter', 'latex');


turbine_center = -9;
Sz = 3;
buffer = 1;
turbineLinewidth = 1.5;
nacelle_length = 0.3;


hold on
for block = 1:3
    contourf(cleaned(block).X, cleaned(block).Y, cleaned(block).u / 8, 100, 'linestyle', 'none')
    clim([0.35, 1])


    % Draw recangles around control volumes
    x = cleaned(block).X(1,:);
    rectangle('position', [min(x), turbine_center - Sz/2, range(x), Sz], 'linewidth', 1)

end

% Draw turbine
% Rotor
plot([0 0], [turbine_center - 0.5, turbine_center + 0.5], 'color', 'black', 'linewidth', turbineLinewidth)
% Nacelle
plot([0, nacelle_length], [turbine_center, turbine_center], 'color', 'black', 'linewidth', turbineLinewidth)



colormap(slanCM('gothic'))
axis equal
set(gca, 'YDir', 'reverse')
xlabel('$x \mathbin{/} D$', 'interpreter', 'latex', 'Fontsize', labelFontSize);
ylabel('$z \mathbin{/} D$', 'interpreter', 'latex', 'Fontsize', labelFontSize);

% Digitized half-width fit
halfwidth_data = readmatrix('/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MKE/OndrejDigitized_HaldWidthFit.csv');

% Reorder
halfwidth_x = halfwidth_data(:,1);
halfwidth_y = halfwidth_data(:,2);
[halfwidth_x_sorted, I] = sort(halfwidth_x, 'ascend');
halfwidth_y_sorted = halfwidth_y(I);


P = plot(halfwidth_x_sorted, halfwidth_y_sorted, 'color', 'black', ...
         'linewidth', 1, 'linestyle', '--');
P.Color(4) = 0.25;


hold off
yticklabels([])
% box on
xlim([-1, 51])
ylim([turbine_center - Sz/2 - buffer, turbine_center + Sz/2 + buffer])
yticks(-12:3:3)



% Save figure
clc; 
save_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/MKE';
figure_name = 'SingleFarmPaper_MKE_ContourCartoon.pdf';
pause(3)
fprintf('Saving figure...\n')
exportgraphics(fig, fullfile(save_folder, figure_name), 'resolution', 600)
clc; fprintf('Figure Saved!\n')













%% Streamwise averaged profiles, Test

block = 1;
Y = cleaned(block).Y;
y = Y(:,1);

 % Average over cropped regions
avg_advection_x = mean(mean_advection(block).x, 2, 'omitnan');
avg_advection_y = mean(mean_advection(block).y, 2, 'omitnan');

avg_transport_x = mean(turbulent_transport(block).x, 2, 'omitnan');
avg_transport_y = mean(turbulent_transport(block).y, 2, 'omitnan');

avg_production_x = mean(production(block).x, 2, 'omitnan');
avg_production_y = mean(production(block).y, 2, 'omitnan');


% Plot
lw = 2;
clc; close all
figure('color', 'white')
hold on

% Advection
plot(avg_advection_x, y, 'linewidth', lw, 'Displayname', 'Advection x')
plot(avg_advection_y, y, 'linewidth', lw, 'Displayname', 'Advection y')
% plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Transport
plot(avg_transport_x, y, 'linewidth', lw, 'Displayname', 'Transport x')
plot(avg_transport_y, y, 'linewidth', lw, 'Displayname', 'Transport y')
% plot(nan, nan, 'color', 'white', 'displayname', ' ')

% Production
plot(avg_production_x, y, 'linewidth', lw, 'Displayname', 'Production x')
plot(avg_production_y, y, 'linewidth', lw, 'Displayname', 'Production y')
% plot(nan, nan, 'color', 'white', 'displayname', ' ')

hold off
set(gca, 'YDir', 'reverse')

legend('box', 'off')
% xlim([-30, 50])
xlim([-130, 130])
xlabel('Advection / Transport / Production')
ylabel('$z / D$', 'interpreter', 'latex')