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

z = 0; % 0 1.5 7.5 9 10.5 
x = 1;  % 1 3 6 8 22.5 48.5


datapath = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/cropped_intantaneous_Natalie';
plane_name = ['xz=' num2str(x) '_' num2str(z)];

bn = 1;
rn = 2; % this is arbitrary for getting the flat array position

[c] = xy2flat(x,z,datapath,bn,rn);

% ----------------------- VELOCITY SIGNAL FIGURE -------------------------
% rn = [1, 2, 3];
[ut] = getvelsignal(x,z,bn,rn,datapath);
ut_mean = mean(ut,2,"omitnan");

rho_step = 10;  % taylor hypothesis time/physical step thing
D = 0.08; % m
fps = 1385;

% Time step in indicies (not time)
t_step = round(((rho_step*D)/ut_mean(1))*fps);

% Time indicies to plot
% t_plot = 2000:2000+t_step;
t_plot = 1:1+t_step;

% Time to plot (units in seconds)
t_step_seconds = t_step / fps;
t_plot_seconds = t_plot / fps;


fig1 = figure('color', 'white');
hold on
grid on
for uu = 1:1
    plot(t_plot_seconds, ut(uu,t_plot(1):t_plot(end)), "Color", lineColor1(3,:), 'LineWidth', 2)
end

xlim([t_plot(1)/fps t_plot(end)/fps])
xlabel('$time$ [s]', 'Interpreter','latex')
ylabel('$u$ [m/s]','Interpreter','latex')

%% ---------------- check the text position of the coordinate on figure

text(1.7,1,['(x,z) = (' num2str(x) ',' num2str(z) ')'])
fontsize(16,"points")
fontname(fig1,"Times New Roman")

%% ---------------- save figure --------------------------------------
% MATLAB figure
% savefig( ...
%     fig1, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_rec.fig']));
% 
% % PNG
% exportgraphics( ...
%     fig1, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_rec.png']), ...
%     'Resolution', 600);