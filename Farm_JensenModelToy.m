%% Jensen / PARK top-hat wake — top-down (x-y) plane
clear; close all; clc

%% ---- Inputs you'll tinker with ----
D  = 120;     % rotor diameter [m]
A0 = 0.25;    % deficit amplitude at x=0 (often 1 - sqrt(1-Ct)) [-]
k  = 0.05;    % wake expansion / recovery rate [-]
% Uinf = NaN;   % set to NaN for normalized U/Uinf, or set e.g. 8 for m/s
Uinf = 7.5;   % set to NaN for normalized U/Uinf, or set e.g. 8 for m/s

x0 = 0;       % wake starts at x0 [m] (often 0 at rotor plane)

% Domain
xMax = 15*D;  % downstream extent [m]
yHalf = 3*D;  % half-width of top-down slice [m]
Nx = 600;     % resolution downstream
Ny = 300;     % resolution lateral

%% ---- Derived quantities ----
R = D/2;

% Grid (top-down plane)
x = linspace(-0.5*D, xMax, Nx);     % include some upstream for context
y = linspace(-yHalf, yHalf, Ny);
[X,Y] = meshgrid(x,y);

% Downstream distance relative to wake start
Xw = X - x0;

% Wake radius expansion (only for x >= x0)
Rw = R + k*Xw;                       % [m]
Rw(Xw < 0) = NaN;                    % undefined upstream

% Jensen top-hat centerline deficit as function of x (only downstream)
def_x = A0 ./ (1 + 2*k*Xw/D).^2;      % DeltaU/Uinf
def_x(Xw < 0) = 0;

% Apply top-hat: deficit is uniform inside |y| <= Rw(x)
inside = (abs(Y) <= Rw) & (Xw >= 0);
defField = def_x .* inside;          % still normalized deficit field

% Convert to velocity field
if isnan(Uinf)
    U = 1 - defField;                % normalized U/Uinf
    cbarLabel = 'U/U_\infty';
else
    U = Uinf * (1 - defField);       % m/s
    cbarLabel = 'U (m/s)';
end

%% ---- Plot ----
figure('Color','w');
% imagesc(x, y, U);
contourf(x/D, y/D, U, 10, 'linestyle', 'none')
set(gca,'YDir','normal');
axis tight
xlabel('x (m)'); ylabel('y (m)');
title('Jensen / PARK top-hat wake (top-down x–y plane)')
cb = colorbar; ylabel(cb, cbarLabel);

hold on
% Wake edge y = ±Rw(x) for x >= x0
xEdge = linspace(x0, xMax, 400);
RwEdge = R + k*(xEdge - x0);
plot(xEdge/D,  RwEdge/D, 'k-', 'LineWidth', 1.5);
plot(xEdge/D, -RwEdge/D, 'k-', 'LineWidth', 1.5);

% Rotor disk projection at x = x0 (just a line segment in top-down view)
plot([x0 x0]/D, [-R R]/D, 'k-', 'LineWidth', 3);

% Helpful annotation
def_end = A0 / (1 + 2*k*(xMax-x0)/D)^2;
text(x0 + 0.02*xMax, 0.90*yHalf, ...
    sprintf('D=%.1f m, A0=%.3f, k=%.3f, deficit@xMax=%.3f', D, A0, k, def_end), ...
    'Color','k','FontWeight','bold');

hold off

%% ---- Quick "physics knobs" you might want to try ----
% 1) Change A0 to represent Ct: A0 = 1 - sqrt(1 - Ct);
% 2) Change k to see wider vs narrower wake growth.
% 3) Change x0 to shift wake start downstream (e.g., near-wake excluded).
