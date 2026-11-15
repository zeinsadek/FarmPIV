clc; close all; clear

% Paths
addpath("/Users/zeinsadek/Documents/MATLAB/MatlabFunctions")
addpath("/Users/zeinsadek/Documents/MATLAB/colormaps")
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps/slanCM')
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/Natalie_functions')

% choose the color map
map = slanCM('bubblegum');
[numColors, bgr] = size(map);
li_color = 4;

% get number of colors
for cc = 1:li_color
    lineColor1(cc,:) = map(round((numColors/li_color)*cc),:);
end

map = slanCM('guppy');
[numColors, bgr] = size(map);
li_color = 4;

% get number of colors
for cc = 1:li_color
    lineColor2(cc,:) = map(round((numColors/li_color)*cc),:);
end
%% choose a coordinate 
% (1-24192) which is the plane dimensions but in array form

datapath = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/cropped_intantaneous_Natalie';

%%% Center turbine is ar x = 0, edge turbine is at x = 9
% Sign is different from figures in the paper but thats okay

% xs = [1, 3, 6, 8, 22.5, 48.5];
xs = [1, 22.5, 48.5];

% center_zs = [0, 3, 9];
% aisle_zs = [1.5, 7.5, 10.5];
% center_labels = {'$T_{IV}$', '$T_{III}$', '$T_{I}$'};
% aisle_labels = {'$T_{III - IV}$', '$T_{I - II}$', 'Aisle'};

% Edge turbine + aisle on either side
center_zs = [9];
aisle_zs = [7.5, 10.5];

center_labels = { '$T_{I}$'};
aisle_labels = {'$T_{I - II}$', 'Aisle'};

% Center turbine plus either side
% center_zs = [0];
% aisle_zs = [1.5, -1.5];
% 
% center_labels = { '$T_{IV}$'};
% aisle_labels = {'$T_{III - IV}$', '$T_{IV - V}$'};
 
center_colors = {'#f5a201', '#ffba42', '#ffd35b'};
aisle_colors = {'#013c58', '#00537a', '#a8e8f9'};


rho_step = 10;  % taylor hypothesis time/physical step thing

u_inf = 8;
D = 0.08;
fps = 1385;

rn = 3;

tickFontSize = 8;
labelFontSize = 10;

lw = 1;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 8, 8]);
tile = tiledlayout(length(xs), 1, 'padding', 'tight', 'tilespacing', 'tight');

for j = 1:length(xs)

    x = xs(j);
    h(j) = nexttile;
    hold on
    ax = gca;
    set(ax, 'FontSize', tickFontSize)
    set(ax, 'TickLabelInterpreter', 'latex');

    % Get plane based on where x location is 
    % BLock 1
    if x >= 1 && x < 4
        bn = 1;
    elseif x >= 4 && x < 7
        bn = 2;
    elseif x >= 7 && x < 10
        bn = 3;
    % Block 2
    elseif x >= 21 && x < 24
        bn = 4;
    elseif x >= 24 && x < 27
        bn = 5;
    elseif x >= 27 && x < 30
        bn = 6;
    % Block 3
    elseif x >= 41 && x < 44
        bn = 7;
    elseif x >= 44 && x < 47
        bn = 8;
    elseif x >= 47 && x < 50
        bn = 9;
    end

    %%% Loop through centers and aisle
    for a = 1:2

        % Centers
        if a == 1
            zs = center_zs;
            colors = center_colors;
            fprintf('Centers\n')
        elseif a == 2
            zs = aisle_zs;
            colors = aisle_colors;
            fprintf('Aisles\n')
        end

        fprintf('Getting signals at x/D = %2.0f from Plane %1.1f\n\n', x, bn)
        for i = 1:length(zs)
            z = zs(i);
        
            % [c] = xy2flat(x,z,datapath,bn,rn);
            fprintf('Getting signal at z/D = %2.1f\n', z)
            
            % ----------------------- VELOCITY SIGNAL FIGURE -------------------------
            % rn = [1, 2, 3];
            [ut] = getvelsignal(x,z,bn,rn,datapath);
            ut_mean = mean(ut, 2, "omitnan");
        
            % Advection time to normalize
            t_advection = D / ut_mean(1);
        
            % Time vector
            time_seconds = 1:size(ut, 2);
            time_seconds = time_seconds - 1;
            time_seconds = time_seconds ./ fps;
        
            % Normalize time by flow cycles
            time_advection = time_seconds ./ t_advection;
    
            % plot(time_advection, fillmissing(ut / u_inf, 'linear'), 'LineWidth', lw, 'color', colors(i,:))
            plot(time_advection, fillmissing(ut / u_inf, 'linear'), 'LineWidth', lw, 'color', colors{i}, 'HandleVisibility', 'off')
            fprintf('\n')
        end
    end

    text(0.25, 1.75, sprintf('$x \\mathbin{/} D = %1.1f$', x), ...
        'interpreter', 'latex', 'fontsize', tickFontSize, ...
        'HorizontalAlignment', 'left')

    hold off
    fprintf('\n\n')

    if j ~= length(xs)
        ax = gca;
        ax.XAxis.Visible = 'off';
    end

    if j == length(xs)
        %%% Legend for colors and centers vs aisles
        % Centers legends
        hold on
        for iii = 1:length(center_zs)
            label = center_labels{iii};
            plot(nan, nan, 'linestyle', '-', 'linewidth', 1, 'Color', center_colors{iii}, 'displayname', label)
        end
        plot(nan, nan, 'color', 'white', 'displayname', ' ')
    
        % Aisle legends
        hold on
        for iii = 1:length(aisle_zs)
            label = aisle_labels{iii};
            plot(nan, nan, 'linestyle', '-', 'linewidth', 1, 'Color', aisle_colors{iii}, 'displayname', label)
        end

        leg = legend('fontsize', tickFontSize, 'interpreter', 'latex', 'Orientation', 'vertical', 'box', 'off');
        leg.Layout.Tile = 'east';
        leg.IconColumnWidth = 10;
    end
end
    

linkaxes(h, 'xy')
xlim([0, rho_step])
ylim([0, 2.0])
xlabel(tile, '$t \cdot \overline{u} \mathbin{/} D$', 'Interpreter','latex', 'fontsize', labelFontSize)
ylabel(tile, '$u \mathbin{/} u_{\infty}$','Interpreter','latex', 'fontsize', labelFontSize)



%% Experimenting with highlighting intermittency

x = 1;
z = 9;

[ut] = getvelsignal(x,z,1,rn,datapath);
ut = fillmissing(ut, 'linear');
% ut = ut(:);
ut_mean = mean(ut, 2, "omitnan");

% Advection time to normalize
t_advection = D / ut_mean(1);

% Time vector
time_seconds = 1:size(ut, 2);
time_seconds = time_seconds - 1;
time_seconds = time_seconds ./ fps;

% Normalize time by flow cycles
time_normalized = time_seconds ./ t_advection;

u_fluc = ut - ut_mean;
u_fluc_normalized = u_fluc / u_inf;


% Gradient of both normalized quantities
u_gradient_normalzied = gradient(u_fluc_normalized, time_normalized);

close all
figure('color', 'white')
hold on
plot(time_normalized, u_gradient_normalzied)
plot(time_normalized, abs(u_gradient_normalzied))
hold off

xlim([0, 10])

% Gradient threshold for intermitency
thresh = 2;
intermittent_ut = ut / u_inf;
intermittent_ut(abs(u_gradient_normalzied) < thresh) = nan;

figure('color', 'white')
hold on
plot(time_normalized, ut / u_inf, 'color', 'black')
P = plot(time_normalized, intermittent_ut, 'color', 'red');
uistack(P, 'top')

xlim([0, 10])
ylim([0, 1])


%% Figure with highlighted intermittency


% xs = [1, 22.5, 48.5];
xs = [3.1];

% Edge turbine + aisle on either side
center_zs = [9];
% aisle_zs = [7.5, 10.5];
aisle_zs = [7.5];

center_labels = { '$T_{I}$'};
aisle_labels = {'$T_{I - II}$', 'Aisle'};

rho_step = 10;

u_inf = 8;
D = 0.08;
fps = 1385;

rn = 1;

tickFontSize = 8;
labelFontSize = 10;

lw = 1;

thresh = 2;

clc; close all
fig = figure('color', 'white', 'units', 'centimeters', 'position', [10, 10, 7.5, 5]);
tile = tiledlayout(length(xs), 1, 'padding', 'tight', 'tilespacing', 'tight');

for j = 1:length(xs)

    x = xs(j);
    h(j) = nexttile;
    hold on
    ax = gca;
    set(ax, 'FontSize', tickFontSize)
    set(ax, 'TickLabelInterpreter', 'latex');

    % Get plane based on where x location is 
    % BLock 1
    if x >= 1 && x < 4
        bn = 1;
    elseif x >= 4 && x < 7
        bn = 2;
    elseif x >= 7 && x < 10
        bn = 3;
    % Block 2
    elseif x >= 21 && x < 24
        bn = 4;
    elseif x >= 24 && x < 27
        bn = 5;
    elseif x >= 27 && x < 30
        bn = 6;
    % Block 3
    elseif x >= 41 && x < 44
        bn = 7;
    elseif x >= 44 && x < 47
        bn = 8;
    elseif x >= 47 && x < 50
        bn = 9;
    end

    %%% Loop through centers and aisle
    for a = 1:2

        % Centers
        if a == 1
            zs = center_zs;
            colors = center_colors;
            fprintf('Centers\n')
        elseif a == 2
            zs = aisle_zs;
            colors = aisle_colors;
            fprintf('Aisles\n')
        end

        fprintf('Getting signals at x/D = %2.0f from Plane %1.1f\n\n', x, bn)
        for i = 1:length(zs)
            z = zs(i);
        
            % [c] = xy2flat(x,z,datapath,bn,rn);
            fprintf('Getting signal at z/D = %2.1f\n', z)
            
            % ----------------------- VELOCITY SIGNAL FIGURE -------------------------
            % rn = [1, 2, 3];
            [ut] = getvelsignal(x,z,bn,rn,datapath);
            ut = fillmissing(ut, 'linear');
            ut_mean = mean(ut, 2, "omitnan");
        
            % Advection time to normalize
            t_advection = D / ut_mean(1);
        
            % Time vector
            time_seconds = 1:size(ut, 2);
            time_seconds = time_seconds - 1;
            time_seconds = time_seconds ./ fps;
        
            % Normalize time by flow cycles
            time_advection = time_seconds ./ t_advection;

            % Get fluctuations
            u_fluc = ut - ut_mean;
            u_fluc_normalized = u_fluc / u_inf;
            
            % Gradient of both normalized quantities
            u_gradient_normalzied = gradient(u_fluc_normalized, time_normalized);

            % Gradient threshold for intermitency
            intermittent_ut = ut / u_inf;
            intermittent_ut(abs(u_gradient_normalzied) < thresh) = nan;
    
            plot(time_advection, ut / u_inf, 'LineWidth', lw, 'color', colors{i}, 'HandleVisibility', 'off')
            % plot(time_advection, intermittent_ut, 'linewidth', lw, 'color', 'red', 'HandleVisibility', 'off')
            fprintf('\n')
        end
    end

    % text(0.25, 1.75, sprintf('$x \\mathbin{/} D = %1.1f$', x), ...
    %     'interpreter', 'latex', 'fontsize', tickFontSize, ...
    %     'HorizontalAlignment', 'left')

    hold off
    fprintf('\n\n')

    if j ~= length(xs)
        ax = gca;
        ax.XAxis.Visible = 'off';
    end

    if j == length(xs)
        %%% Legend for colors and centers vs aisles
        % Centers legends
        hold on
        for iii = 1:length(center_zs)
            label = center_labels{iii};
            plot(nan, nan, 'linestyle', '-', 'linewidth', 1, 'Color', center_colors{iii}, 'displayname', label)
        end
        % plot(nan, nan, 'color', 'white', 'displayname', ' ')
    
        % Aisle legends
        hold on
        for iii = 1:length(aisle_zs)
            label = aisle_labels{iii};
            plot(nan, nan, 'linestyle', '-', 'linewidth', 1, 'Color', aisle_colors{iii}, 'displayname', label)
        end

        leg = legend('fontsize', tickFontSize, 'interpreter', 'latex', 'Orientation', 'horizontal', 'box', 'off');
        leg.Layout.Tile = 'north';
        leg.IconColumnWidth = 10;
    end
end
    

linkaxes(h, 'xy')
xlim([0, rho_step])
ylim([0.2, 1.2])
yticks(0.2:0.4:1.2)
xlabel(tile, '$t \cdot \overline{u} \mathbin{/} D$', 'Interpreter','latex', 'fontsize', labelFontSize)
ylabel(tile, '$u \mathbin{/} u_{\infty}$','Interpreter','latex', 'fontsize', labelFontSize)

% Save figure
figure_folder = '/Users/zeinsadek/Desktop/Experiments/Farm/SingleFarmPaper/VelocitySignal';
filename = 'SingleFarm_VelocitySignal.pdf';
pause(3)
exportgraphics(fig, fullfile(figure_folder, filename), 'resolution', 600)
close all







