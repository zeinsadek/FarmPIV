clear
% path where all the matlab functions needed are
addpath(genpath('/Volumes/ATHENA/MATLAB_CODE'));

%% things to change

% coordinates
z = 10.5;        % 0 1.5 7.5 9 10.5 
x = 48.5;          % 1 3 6 8 22.5 48.5

% block for that coordinate. If this is a mismatch, it will return an error
bn = 9;

% path for data and where to save the results
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';
pathout = '/Volumes/ATHENA/FARM/POSTPROC/blocks/SingleFarm/castaing';

% colormap from the slanCM folder
map = slanCM('bubblegum');

%% ------------------ PREAMBLE ---------------------------------------
% rho values that we want (this is all we have rn so don't go crazy)
rho_input = [0.08, 0.1, 0.5, 0.75, 1, 2, 3, 4, 5, 7, 10, 18, 63];

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
fig3 = figure()
axes('XScale', 'linear', 'YScale', 'log')
hold on
grid on
for i = 1:length(rho_input)
    plot(data(i).pdf(c).xpdf, (data(i).pdf(c).pdf)*10^i, '-o','MarkerEdgeColor',lineColor1(i,:),'MarkerFaceColor',lineColor1(i,:),'MarkerSize',3,'Color',lineColor1(i,:))
end
% xlim([-1.5 1.5])
ylabel('$PDF$','Interpreter','latex')
xlabel('$\delta u_\tau$', 'Interpreter', 'latex')
xlim([-6, 6])
text(3.5,10^14,['(x,z) = (' num2str(x) ',' num2str(z) ')'])
h=legend('0.08','0.1','0.5','0.75','1','2','3','4','5','7','10','18','63','Location','eastoutside')
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

beep

