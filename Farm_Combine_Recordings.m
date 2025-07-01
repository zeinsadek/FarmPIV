% Code to combine all three recordings of the TR-PIV 
% to compute the means and stresses

% Zein Sadek
% Portland State University

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/Farm_Functions');
addpath('/Users/zeinsadek/Desktop/Experiments/PIV/Processing/readimx-v2.1.8-osx');
addpath('/Users/zeinsadek/Documents/MATLAB/colormaps');
fprintf("All Paths Imported...\n\n");

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DATA + SAVE PATHS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

experiment = 'Farm2Farm_40D_Gap';
plane = "Plane_9";

results_path = fullfile("/Users/zeinsadek/Library/Mobile Documents/com~apple~CloudDocs/Data/Farm");
data_path = fullfile(results_path, "data", experiment);
combined_path = fullfile(results_path, "combined", experiment);
combined_name = strcat(plane, "_COMBINED.mat");

if ~exist(combined_path, "dir")
    mkdir(combined_path)
end

% Find how many recordings are available in directory
files = dir(data_path);
names = {files.name};
names = erase(names, '.mat');
names = names(~contains(names, {'._', '..', '.'}));
recordings = sum(contains(names, plane));

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for r = 1:recordings
    fprintf("Importing Recording %1.0f\n\n", r);
    name = strcat(plane, "_Recording_", num2str(r), "_DATA.mat");
    tmp = load(fullfile(data_path, name));
    data(r) = tmp.output;
end

clear files names
clear tmp name r
clc;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMBINE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

combined.U = cat(3, data(:).U);
combined.V = cat(3, data(:).V);
combined.X = data(1).X;
combined.Y = data(1).Y;
combined.D = length(combined.V);
combined.R = recordings;

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPUTE MEANS AND STRESSES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_file = fullfile(combined_path, combined_name);
combined_means = data2means(save_file, combined);