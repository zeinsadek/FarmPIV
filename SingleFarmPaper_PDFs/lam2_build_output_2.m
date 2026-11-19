%% -------------------- THINGS TO CHANGE -----------------------------

clc; close all; clear

% One row per requested [x/D, z/D] pair, in plotting order.
z = 10.5;
% xz_pairs = [
%      1, z
%      2, z
%      3, z
%      4, z
%      5, z
%      6, z
%      7, z
%      8, z
%      9, z
% 
%      21, z
%      22, z
%      23, z
%      24, z
%      25, z
%      26, z
%      27, z
%      28, z
%      29, z
% 
%      41, z
%      42, z
%      43, z
%      44, z
%      45, z
%      46, z
%      47, z
%      48, z
%      49, z
% ];

xz_pairs = [
     2.5, z
     5.5, z
     8.5, z

     22.5, z
     25.5, z
     28.5, z

     42.5, z
     45.5, z
     48.5, z
];

rho_input = [0.08, 0.1, 0.5, 0.75, 1, 2, 3, 4, 5, 7, 10, 18, 63];
fit_range = [0.1, 3];   % Inclusive bounds in rho/D
ver = 8;
rn = 1;

pathin = '/Volumes/ATHENA/FARM/POSTPROC/blocks/SingleFarm/castaing';
datapath = '/Volumes/ATHENA/FARM/DATA/blocks/SingleFarm';
output_file = 'lam2_output.mat';  % Saved in the current MATLAB folder

%% -------------------- FIND THE PLANE AND GRID POINT -----------------

validateattributes(xz_pairs, {'numeric'}, ...
    {'real', 'finite', '2d', 'nonempty', 'ncols', 2});
validateattributes(rho_input, {'numeric'}, ...
    {'real', 'finite', 'positive', 'vector', 'increasing'});
validateattributes(fit_range, {'numeric'}, ...
    {'real', 'finite', 'positive', 'vector', 'numel', 2, 'increasing'});

plane_start = [1, 4, 7, 21, 24, 27, 41, 44, 47];
plane_width = 3;
n_cases = size(xz_pairs, 1);
plane = zeros(n_cases, 1);
flat_index = zeros(n_cases, 1);

% Resolve all planes before accessing the data.
for j = 1:n_cases
    x = xz_pairs(j, 1);
    bn = find(x >= plane_start & x < plane_start + plane_width, 1);
    if isempty(bn)
        error('x/D = %g falls outside the measured plane intervals.', x);
    end
    plane(j) = bn;
end

rho_input = rho_input(:).';
fit_mask = rho_input >= fit_range(1) & rho_input <= fit_range(2);
if nnz(fit_mask) < 2
    error('The fit range must contain at least two sampled rho/D values.');
end

for j = 1:n_cases
    flat_index(j) = xy2flat(xz_pairs(j, 1), xz_pairs(j, 2), ...
        datapath, plane(j), rn);
end

%% -------------------- LOAD LAMBDA-SQUARED DATA ----------------------

lambda_sq = nan(n_cases, numel(rho_input));
planes_used = unique(plane, 'stable').';

for bn = planes_used
    case_idx = find(plane == bn);
    point_idx = flat_index(case_idx);

    for i = 1:numel(rho_input)
        rho = rho_input(i);
        file_name = fullfile(pathin, ['rho_' num2str(rho, 3)], ...
            ['v' num2str(ver)], ...
            sprintf('Plane_%d_Recording_all_rho_%g.mat', bn, rho));

        fprintf('Loading %s\n', file_name);
        data = load(file_name, 'lambda_sq');
        if ~isfield(data, 'lambda_sq')
            error('Variable lambda_sq was not found in %s.', file_name);
        end

        % MATLAB linear indexing matches the original column-major reshape.
        values = data.lambda_sq(point_idx);
        lambda_sq(case_idx, i) = values(:);
    end
end

%% -------------------- FIT LINES AND SAVE OUTPUT ---------------------

output = struct();
output.rho = rho_input;
output.fit_range = fit_range;
output.fit_mask = fit_mask;
output.version = ver;
output.recording = rn;

% Save the original plot settings so the plotting script only needs the MAT file.
map = slanCM('bubblegum', max(256, n_cases));
color_idx = round((1:n_cases) * size(map, 1) / n_cases);
output.plot.x_scale = 'log';
output.plot.x_label = '$\rho/D$';
output.plot.y_label = '$\lambda^2$';
output.plot.y_lim = [-0.05, 0.4];
output.plot.n_cols = min(3, n_cases);
output.plot.fit_color = [1, 0, 0];
output.plot.marker_size = 5;
output.plot.font_name = 'Times New Roman';
output.plot.font_size = 14;
output.plot.label_xy = [1, 0.30; 1, 0.25];

for j = 1:n_cases
    y = lambda_sq(j, :);
    valid_fit = fit_mask & isfinite(y);
    rho_fit = rho_input(valid_fit);
    log_rho_fit = log10(rho_fit);

    p = [NaN, NaN];
    y_fit = [];
    if nnz(valid_fit) >= 2
        p = polyfit(log_rho_fit, y(valid_fit), 1);
        y_fit = polyval(p, log_rho_fit);
    else
        warning('Fewer than two finite fit samples at (x,z) = (%g,%g).', ...
            xz_pairs(j, 1), xz_pairs(j, 2));
    end

    output.cases(j).x = xz_pairs(j, 1);
    output.cases(j).z = xz_pairs(j, 2);
    output.cases(j).plane = plane(j);
    output.cases(j).flat_index = flat_index(j);
    output.cases(j).lambda_sq = y;
    output.cases(j).fit_mask = valid_fit;
    output.cases(j).fit.coefficients = p;  % [slope, intercept]
    output.cases(j).fit.slope = p(1);  % d(lambda^2)/d(log10(rho/D))
    output.cases(j).fit.rho = rho_fit;
    output.cases(j).fit.log10_rho = log_rho_fit;
    output.cases(j).fit.lambda_sq_data = y(valid_fit);
    output.cases(j).fit.lambda_sq_line = y_fit;
    output.cases(j).color = map(color_idx(j), :);
end


save_folder = '/Users/zeinsadek/Desktop/Experiments/PIV/Processing/Farm/SingleFarmPaper_PDFs/CastaingMatfiles';
filename = sprintf('CastaingPlotting_z_%1.0f.mat', z);
save(fullfile(save_folder, filename), 'output');

