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


%% Save signal into matfile

datapath = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/cropped_intantaneous_Natalie';

%%% Center turbine is ar z = 0, edge turbine is at z = 9
% Sign is different from figures in the paper but thats okay
x = 48.5;
z = 0;

u_inf = 8;
D = 0.08;
fps = 1385;

rn = 1;

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

output.time_seconds = time_seconds;
output.time_normalized = time_advection;
output.advection_time = t_advection;
output.u_mean = ut_mean;

output.signal = ut;
output.signal_normalized = ut / u_inf;

output.x = x;
output.z = z;
output.plane = bn;


% Save figure
% save_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_VelocitySignal/SingleSignalMatfiles';
% filename = sprintf('SingleFarm_VelocitySignal_x_%1.0f_z_%1.0f.mat', x, z);
% save(fullfile(save_folder, filename), 'output');







