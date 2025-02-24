%%% Check Convergence Across Multiple Time-Resolved Recordings
% Ondenburg, DE
% Zein Sadek

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/colormaps');
fprintf("All Paths Imported...\n\n");
% addpath('G:/Other computers/Zein MacBook Pro/Farm/Farm_Functions/');
% addpath('C:/Users/Zein/Documents/MATLAB/readimx-v2.1.8-win64/');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERGENCE CALCULATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;
data_folder    = '/Volumes/Zein_PIV_4/Oldenburg_Convergence/data';
num_recordings = 3;
PIV_plane      = 1;
step           = 10;

r = 20;
c_center_wake = 403;
c_aisle = 307;
c_edge = 130;

for i = 1:num_recordings
case_name  = strcat('Plane_', num2str(PIV_plane), '_Recording_', num2str(i));
fprintf(strcat(case_name, '\n'));

data_path = fullfile(data_folder, strcat(case_name, '_DATA.mat'));
data      = load(data_path);
data      = data.output;

N          = data.D;
groups     = 1:step:N;
avgs_wake  = zeros(size(groups));
avgs_aisle = zeros(size(groups));
avgs_edge  = zeros(size(groups));
avgs_all   = zeros(size(groups));

c = 1;
for j = step:step:N

    % Print Progress.
    progressbarText(j/N);

    % Wake
    avgs_wake(c) = mean(data.U(r, c_center_wake, 1:j), 'omitnan');

    % Aisle
    avgs_aisle(c) = mean(data.U(r, c_aisle, 1:j), 'omitnan');

    % Edge
    avgs_edge(c) = mean(data.U(r, c_edge, 1:j), 'omitnan');

    % All
    avgs_all(c) = mean(data.U(:, :, 1:j), 'all', 'omitnan');
    c = c + 1;
end

convergence.(case_name).wake  = avgs_wake;
convergence.(case_name).aisle = avgs_aisle;
convergence.(case_name).edge  = avgs_edge;
convergence.groups = groups;

end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AVERAGING ACROSS MULTIPLE RECORDINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data sets
rec1 = load('/Volumes/Zein_PIV_4/Oldenburg_Convergence/data/Plane_1_Recording_1_DATA.mat');
rec1 = rec1.output;
rec1_u = rec1.U;

rec2 = load('/Volumes/Zein_PIV_4/Oldenburg_Convergence/data/Plane_1_Recording_2_DATA.mat');
rec2 = rec2.output;
rec2_u = rec2.U;

rec3 = load('/Volumes/Zein_PIV_4/Oldenburg_Convergence/data/Plane_1_Recording_3_DATA.mat');
rec3 = rec3.output;
rec3_u = rec3.U;

% Concatenate
big_data = cat(3, rec1_u, rec2_u, rec3_u);

%%

% Compute Convergence
step       = 100;
D          = size(big_data, 3);
groups     = 1:step:D;
avgs_wake  = zeros(size(groups));
avgs_aisle = zeros(size(groups));
avgs_edge  = zeros(size(groups));
avgs_all   = zeros(size(groups));

c = 1;
for j = step:step:D

    % Print Progress.
    progressbarText(j/D);

    % Wake
    avgs_wake(c) = mean(big_data(r, c_center_wake, 1:j), 'omitnan');

    % Aisle
    avgs_aisle(c) = mean(big_data(r, c_aisle, 1:j), 'omitnan');

    % Edge
    avgs_edge(c) = mean(big_data(r, c_edge, 1:j), 'omitnan');

    % All
    avgs_all(c) = mean(big_data(:, :, 1:j), 'all', 'omitnan');
    c = c + 1;
end

convergence.all.wake   = avgs_wake;
convergence.all.aisle  = avgs_aisle;
convergence.all.edge   = avgs_edge;
convergence.all.groups = groups;




%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERGENCE PLOTTING 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ls = 2;
fig = figure('position', [300, 500, 1000, 400]);

hold on
for i = 1:num_recordings
case_name  = strcat('Plane_', num2str(PIV_plane), '_Recording_', num2str(i));
wake = convergence.(case_name).wake;
aisle = convergence.(case_name).aisle;
edge = convergence.(case_name).edge;

if i == 1
    color = 'red';
elseif i == 2
    color = 'green';
elseif i == 3
    color = 'blue';
end

plot(convergence.groups, wake, 'color', color, 'linestyle', '-', 'linewidth', ls, 'displayname', case_name);
plot(convergence.groups, aisle, 'color', color, 'linestyle', '--', 'linewidth', ls, 'displayname', case_name);
plot(convergence.groups, edge, 'color', color, 'linestyle', '-.', 'linewidth', ls, 'displayname', case_name);

end

% 
% plot(convergence.all.groups, convergence.all.wake, 'linestyle', '-', 'linewidth', ls, 'color', 'black', 'displayname', 'none');
% plot(convergence.all.groups, convergence.all.aisle, 'linestyle', '--', 'linewidth', ls, 'color', 'black', 'displayname', 'none');
% plot(convergence.all.groups, convergence.all.edge, 'linestyle', '-.', 'linewidth', ls, 'color', 'black', 'displayname', 'none');
% 



hold off
title('Convergence')
xlabel('Num Images');
ylabel('Mean [m/s]');
% xlim([step, N - 2 * step]);
ylim([0, 10]);


%%

ax = figure();
hold on

plot(convergence.all.groups, convergence.all.wake, 'linestyle', '-', 'linewidth', ls, 'color', 'black', 'displayname', 'none');
plot(convergence.all.groups, convergence.all.aisle, 'linestyle', '--', 'linewidth', ls, 'color', 'black', 'displayname', 'none');
plot(convergence.all.groups, convergence.all.edge, 'linestyle', '-.', 'linewidth', ls, 'color', 'black', 'displayname', 'none');

xline(N, 'r');
xline(2*N, 'r');

yline(convergence.all.wake(end-10))
yline(convergence.all.aisle(end-10))
yline(convergence.all.edge(end-10))
hold off
xlim([step, D - 2*step])
xlabel('Number of Images')
ylabel('u [m/s]')
title('Convergence: 3 Combined, Independent, Time-Resolved Measurements')























