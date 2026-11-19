clear

%% ------------ THINGS TO CHANGE ---------------------------

% ---------- saving location --------------------------------------
pathout = '/Volumes/ATHENA/FARM/POSTPROC/blocks/SingleFarm/castaing/lam2_vrho';
map = slanCM('bubblegum');


%% -------------------------- RHO VALUES -----------------------------
% (1-24192) which is the plane dimensions but in array form
rho_input = [0.08, 0.1, 0.5, 0.75, 1, 2, 3, 4, 5, 7, 10, 18, 63];
ver = 8;
XX = [1, 3, 6, 8, 22.5, 48.5];
ZZ = [0, 1.5, 7.5, 9, 10.5];

% ------------------- CHOOSE THE COLOR MAP -----------------------------
[numColors, bgr] = size(map);
li_color = length(rho_input);

% get number of colors
for cc = 1:li_color
    lineColor1(cc,:) = map(round((numColors/li_color)*cc),:);
end

map = slanCM('bubblegum');
[numColors, bgr] = size(map);
li_color = 6;

% get number of colors
for cc = 1:li_color
    lineColor2(cc,:) = map(round((numColors/li_color)*cc),:);
end

%% -------------------- RETRIEVE DATA ----------------------------------
% this is for every coordinate. it takes some time to run. 
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 1;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [1]; % xx(1)
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);

for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx1(i) = rho_input(i);
        lam2_y1(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

% ----------------------------- block 1 --------------------------- 
clearvars -except rho_dx1 lam2_y1 lineColor1 lineColor2 rho_input ver
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 1;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [3];
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);

for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx11(i) = rho_input(i);
        lam2_y11(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

% ----------------------- block 2 ----------------------------------

clearvars -except ver rho_dx1 lam2_y1 lineColor1 lineColor2 rho_input rho_dx11 lam2_y11
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 2;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [6];
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);
%

for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx2(i) = rho_input(i);
        lam2_y2(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

% ----------------------- block 3 ----------------------------------
clearvars -except ver rho_dx1 lam2_y1 rho_dx11 lam2_y11 rho_dx2 lam2_y2 lineColor1 lineColor2 rho_input
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 3;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [8];
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);
%

for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx3(i) = rho_input(i);
        lam2_y3(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

% ------------------------------- block 4 -----------------------------
clearvars -except ver rho_dx1 lam2_y1 rho_dx11 lam2_y11 rho_dx2 lam2_y2 rho_dx3 lam2_y3 lineColor1 lineColor2 rho_input
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 4;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [22.5];
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);
%
for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx4(i) = rho_input(i);
        lam2_y4(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

% -------------------- block 9 ----------------------------------------
clearvars -except ver rho_dx1 lam2_y1 rho_dx11 lam2_y11 rho_dx2 lam2_y2 rho_dx3 lam2_y3 rho_dx4 lam2_y4 lineColor1 lineColor2 rho_input
DataType = 'blocks';
FarmCase = 'SingleFarm';

% for one plane 
pn = 9;

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

zz = [0, 1.5, 7.5, 9, 10.5];
xx = [48.5];
bn = pn;
rn = 1;
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';

for j = 1:length(xx)
    for jj = 1:length(zz)
        [c_loc(j,jj)] = xy2flat(xx(j),zz(jj),datapath,bn,rn);
    end
end
cx = size(c_loc,1);
cy = size(c_loc,2);

c_loc_flat = reshape(c_loc,[1, cx*cy]);

% get the lam^2 data
for ii = 1:length(c_loc_flat)
    for i = 1:length(rho_input)
        clear Nx Ny lambda_sq_flat
        Nx = size(data(i).lambda_sq,1);
        Ny = size(data(i).lambda_sq,2);
        lambda_sq_flat = reshape(data(i).lambda_sq,[1, Nx*Ny]);
        rho_dx9(i) = rho_input(i);
        lam2_y9(ii,i) = lambda_sq_flat(c_loc_flat(ii));
    end
end

%% ------------------------ LINE FITS -------------------------------
% now that you have the lam^2 for each x and z coordinate, fit a line to
% the linear part. This is based on the range that you give. Right now i
% have it fitting from rho/D = 0.1:3. 

fit_range = 2:7; % index values

% for lam2_y1
for i = 1:size(lam2_y1,1)
    [PL1(i,:), SL1(i)] = polyfit(log10(rho_dx1(fit_range)),lam2_y1(i,fit_range),1);
    yLinear1(i,:) = polyval(PL1(i,:),log10(rho_dx1(fit_range)));
end 

% for lam2_y11
for i = 1:size(lam2_y11,1)
    [PL11(i,:), SL11(i)] = polyfit(log10(rho_dx11(fit_range)),lam2_y11(i,fit_range),1);
    yLinear11(i,:) = polyval(PL11(i,:),log10(rho_dx11(fit_range)));
end 

% for lam2_y2
for i = 1:size(lam2_y2,1)
    [PL2(i,:), SL2(i)] = polyfit(log10(rho_dx2(fit_range)),lam2_y2(i,fit_range),1);
    yLinear2(i,:) = polyval(PL2(i,:),log10(rho_dx2(fit_range)));
end 

% for lam2_y3
for i = 1:size(lam2_y3,1)
    [PL3(i,:), SL3(i)] = polyfit(log10(rho_dx3(fit_range)),lam2_y3(i,fit_range),1);
    yLinear3(i,:) = polyval(PL3(i,:),log10(rho_dx3(fit_range)));
end 

% for lam2_4
for i = 1:size(lam2_y4,1)
    [PL4(i,:), SL4(i)] = polyfit(log10(rho_dx4(fit_range)),lam2_y4(i,fit_range),1);
    yLinear4(i,:) = polyval(PL4(i,:),log10(rho_dx4(fit_range)));
end 

% for lam2_y9
for i = 1:size(lam2_y9,1)
    [PL9(i,:), SL9(i)] = polyfit(log10(rho_dx9(fit_range)),lam2_y9(i,fit_range),1);
    yLinear9(i,:) = polyval(PL9(i,:),log10(rho_dx9(fit_range)));
end 

% 
% 
% 
% figure()
% scatter(rho_dx3, lam2_y3(1,:),'o')
% hold on
% plot(rho_dx3(fit_range),yLinear,'-')
% grid on
% set(gca,"XScale","log")


%% ----------- FIGURES -------------------------------------------
% this figure is for z = 0, 1.5, 7.5, 9, 10.5 with x = 1, 3, 6, 8, 22.5, 48.5
z = 5; % index value that indicates the array above

fig3 = figure()
t3 = tiledlayout(2,3)
ax1 = nexttile
hold on
plot(rho_dx1, lam2_y1(z,:),'-o','Color',lineColor2(1,:),'MarkerEdgeColor',lineColor2(1,:),'MarkerFaceColor',lineColor2(1,:),'MarkerSize',5)
grid on
plot(rho_dx1(fit_range),yLinear1(z,:), '-','Color',"r")

ax2 = nexttile
hold on
plot(rho_dx11, lam2_y11(z,:),'-o','Color',lineColor2(2,:),'MarkerEdgeColor',lineColor2(2,:),'MarkerFaceColor',lineColor2(2,:),'MarkerSize',5)
grid on
plot(rho_dx11(fit_range),yLinear11(z,:), '-','Color',"r")

ax3 = nexttile
hold on
plot(rho_dx2, lam2_y2(z,:),'-o','Color',lineColor2(3,:),'MarkerEdgeColor',lineColor2(3,:),'MarkerFaceColor',lineColor2(3,:),'MarkerSize',5)
grid on
plot(rho_dx2(fit_range),yLinear2(z,:), '-','Color',"r")

ax4 = nexttile
hold on
plot(rho_dx3, lam2_y3(z,:),'-o','Color',lineColor2(4,:),'MarkerEdgeColor',lineColor2(4,:),'MarkerFaceColor',lineColor2(4,:),'MarkerSize',5)
grid on
plot(rho_dx3(fit_range),yLinear3(z,:), '-','Color',"r")

ax5 = nexttile
hold on
plot(rho_dx4, lam2_y4(z,:),'-o','Color',lineColor2(5,:),'MarkerEdgeColor',lineColor2(5,:),'MarkerFaceColor',lineColor2(5,:),'MarkerSize',5)
grid on
plot(rho_dx4(fit_range),yLinear4(z,:), '-','Color',"r")

ax6 = nexttile
hold on
plot(rho_dx9, lam2_y9(z,:),'-o','Color',lineColor2(6,:),'MarkerEdgeColor',lineColor2(6,:),'MarkerFaceColor',lineColor2(6,:),'MarkerSize',5)
grid on
plot(rho_dx9(fit_range),yLinear9(z,:), '-','Color',"r")

%------------- general figure stuff: axes, title, legend, limits, etc
ax1.XScale = "log";
ax2.XScale = "log";
ax3.XScale = "log";
ax4.XScale = "log";
ax5.XScale = "log";
ax6.XScale = "log";

ax1.YLim = [-0.05 0.4];
ax2.YLim = [-0.05 0.4];
ax3.YLim = [-0.05 0.4];
ax4.YLim = [-0.05 0.4];
ax5.YLim = [-0.05 0.4];
ax6.YLim = [-0.05 0.4];

t3.XLabel.String = '$\rho/D$';
t3.XLabel.Interpreter = 'latex';
t3.YLabel.String = '$\lambda^2$';
t3.YLabel.Interpreter = 'latex';

text(ax1,1,0.3,['(x,z) = (' num2str(XX(1)) ',' num2str(ZZ(z)) ')'])
text(ax2,1,0.3,['(x,z) = (' num2str(XX(2)) ',' num2str(ZZ(z)) ')'])
text(ax3,1,0.3,['(x,z) = (' num2str(XX(3)) ',' num2str(ZZ(z)) ')'])
text(ax4,1,0.3,['(x,z) = (' num2str(XX(4)) ',' num2str(ZZ(z)) ')'])
text(ax5,1,0.3,['(x,z) = (' num2str(XX(5)) ',' num2str(ZZ(z)) ')'])
text(ax6,1,0.3,['(x,z) = (' num2str(XX(6)) ',' num2str(ZZ(z)) ')'])

text(ax1,1,0.25,['m = ' num2str(PL1(z,1))])
text(ax2,1,0.25,['m = ' num2str(PL11(z,1))])
text(ax3,1,0.25,['m = ' num2str(PL2(z,1))])
text(ax4,1,0.25,['m = ' num2str(PL3(z,1))])
text(ax5,1,0.25,['m = ' num2str(PL4(z,1))])
text(ax6,1,0.25,['m = ' num2str(PL9(z,1))])

fontsize(14,"points")
fontname(fig3,"Times New Roman")

% MATLAB figure
savefig( ...
    fig3, ...
    fullfile( ...
    pathout, ...
    ['z' num2str(ZZ(z)) '.fig']));

% PNG
exportgraphics( ...
    fig3, ...
    fullfile( ...
    pathout, ...
    ['z' num2str(ZZ(z)) '.png']), ...
    'Resolution', 600);
























