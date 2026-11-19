clear
% path where all the matlab functions needed are
% addpath(genpath('/Volumes/ATHENA/MATLAB_CODE'));

%% things to change

% coordinates
% z = 10.5;        % 0 1.5 7.5 9 10.5 

%%% Center turbine is ar x = 0, edge turbine is at x = 9
% Sign is different from figures in the paper but thats okay

z = 1.5;        % 0 1.5 7.5 9 10.5 
x = 48.5;          % 1 3 6 8 22.5 48.5

% block for that coordinate. If this is a mismatch, it will return an error
% bn = 3;

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


% path for data and where to save the results
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';
pathout = '/Volumes/ATHENA/FARM/POSTPROC/blocks/SingleFarm/castaing';

% colormap from the slanCM folder
map = slanCM('bubblegum');

%% ------------------ PREAMBLE ---------------------------------------
% rho values that we want (this is all we have rn so don't go crazy)
rho_input = [0.08, 0.1, 0.5, 0.75, 1, 2, 3, 4, 5, 7, 10, 18, 63];
% rho_input = [0.08, 0.1, 0.5];

% name that will be used for saving data
plane_name = ['xz=' num2str(x) '_' num2str(z)];
rn = 1; % this is arbitrary for getting the flat array position

% this function gets the flattened index value of the coordinate you chose
[c] = xy2flat(x,z,datapath,bn,rn);

%% ---------------------- COLOR MAP -------------------------------------
[numColors, bgr] = size(map);
li_color = length(rho_input);

% get number of colors
for cc = 1:li_color
    lineColor1(cc,:) = map(round((numColors/li_color)*cc),:);
end

%% -------------------- RETRIEVE DATA ----------------------------------

DataType = 'blocks';
FarmCase = 'SingleFarm';
ver = 8;

% for one plane 
pn = bn;
% rn = 1;
for i = 1:length(rho_input)
    pathin = ['/Volumes/ATHENA/FARM/POSTPROC/', DataType filesep FarmCase filesep 'castaing/rho_' num2str(rho_input(i),3) filesep 'v' num2str(ver)]; 
    file_search = sprintf( ...
        'Plane_%d_Recording_all_rho_%g.mat', ...
        pn, rho_input(i));

    files = dir(fullfile(pathin, file_search));

    if isempty(files)
        error('No matching data file found.');
    elseif numel(files) > 1
        warning('Multiple matching files found. Using the first file.');
    end

    file_name = fullfile( ...
        files(1).folder, ...
        files(1).name);

    fprintf('\nLoading:\n%s\n\n', file_name);

    data(i) = load(file_name,'pdf','lambda_sq');
end

%% ------------------ GAUSSIAN FITS ------------------------------------
% flatten the lam2
% lx = size(data(1).lambda_sq,1);
% ly = size(data(1).lambda_sq,2);
% lam2flat = reshape(data(1).lambda_sq,[1, lx*ly]);
% 
% % get the pdf at that coordinate 
% xpdf = data(1).pdf(c).xpdf;
% ypdf = data(1).pdf(c).pdf;
% 
% f = fit(xpdf, ypdf, 'gauss1');

%% ------------------- plot the pdf at this coordinate


edge_buffer = 468;

clc; close all
fig3 = figure('color', 'white');
axes('XScale', 'linear', 'YScale', 'log')
hold on
grid on
for i = 1:length(rho_input)

    tmp_x = data(i).pdf(c).xpdf;
    tmp_y = data(i).pdf(c).pdf;

    tmp_y(1:edge_buffer) = nan;
    tmp_y(end - edge_buffer:end) = nan;

    plot(tmp_x, (tmp_y)*10^i, '-o', ...
         'MarkerEdgeColor',lineColor1(i,:), 'MarkerFaceColor',lineColor1(i,:), ...
         'MarkerSize',3,'Color',lineColor1(i,:))
end
% xlim([-1.5 1.5])
ylabel('$PDF$','Interpreter','latex')
xlabel('$\delta u_\tau$', 'Interpreter', 'latex')
xlim([-6, 6])
ylim([0, 1E14])
% text(-5,10^12,['(x,z) = (' num2str(x) ',' num2str(z) ')'], 'fontsize', 4)
title(['(x,z) = (' num2str(x) ',' num2str(z) ')'], 'fontsize', 8)
h=legend('0.08','0.1','0.5','0.75','1','2','3','4','5','7','10','18','63','Location','eastoutside', 'box', 'off')
h.IconColumnWidth = 15;
fontsize(16,"points")
fontname(fig3,"Times New Roman")
h.FontSize = 12;
h.Interpreter = 'latex';



%% ------------------------- SAVING --------------------------------------
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
% 
% % MATLAB figure
% savefig( ...
%     fig2, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_lam.fig']));
% 
% % PNG
% exportgraphics( ...
%     fig2, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_lam.png']), ...
%     'Resolution', 600);

% MATLAB figure
% savefig( ...
%     fig3, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_pdf.fig']));
% 
% % PNG
% exportgraphics( ...
%     fig3, ...
%     fullfile( ...
%     pathout, ...
%     [plane_name '_pdf.png']), ...
%     'Resolution', 600);

% beep
% 1. Read the MP3 file
[y, Fs] = audioread('/Users/zeinsadek/Desktop/Sounds/Fahh.mp3');

% 2. Play the audio using sound (or soundsc for auto-scaling)
sound(y, Fs);
