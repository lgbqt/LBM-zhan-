clear;
clc;
close all;

%% ============================================================
%  Zhan 2023 droplet impact
%  Spreading radius post-processing
%
%  Definition used here:
%
%  At y = h + W,
%  find phi = 0.5 intersections.
%
%  The outermost left/right intersections are used:
%
%        r = (xR - xL)/2
%
%  IMPORTANT:
%  This is OUR reproducible post-processing definition.
%  The paper does not explicitly give its extraction algorithm.
%% ============================================================

files = {
    'droplet_phi_iso_step_00400.dat'
    'droplet_phi_iso_step_00800.dat'
    'droplet_phi_iso_step_01600.dat'
    'droplet_phi_iso_step_03200.dat'
    'droplet_phi_iso_step_06400.dat'
};

T = [
    0.1
    0.2
    0.4
    0.8
    1.6
];

%% ============================================================
% Parameters
%% ============================================================

NX = 1500;
NY = 500;

dx = 1.0;

D = 200.0;

h = 50.0;

W = 5.0;

phiInterface = 0.5;

% Probe height
yProbe = h + W;

%% ============================================================
% Coordinates
%% ============================================================

x = -750 + (0.5 : 1 : 1499.5);

y = 0.5 : 1 : 499.5;

%% ============================================================
% Output arrays
%% ============================================================

rNumerical = zeros(length(T),1);

rOverD = zeros(length(T),1);

rTheoryOverD = zeros(length(T),1);

relativeDifference = zeros(length(T),1);

numberOfCrossings = zeros(length(T),1);

%% ============================================================
% Loop over files
%% ============================================================

for n = 1:length(T)

    phi = readmatrix(files{n});

    %% --------------------------------------------------------
    % Interpolate phi onto y = yProbe
    %% --------------------------------------------------------

    upperIndex = find(y > yProbe,1);

    lowerIndex = upperIndex - 1;

    y1 = y(lowerIndex);
    y2 = y(upperIndex);

    alpha = ...
        (yProbe - y1) ...
        / ...
        (y2 - y1);

    phiLine = ...
        (1.0-alpha) * phi(lowerIndex,:) ...
        + ...
        alpha * phi(upperIndex,:);

    %% --------------------------------------------------------
    % Find phi = 0.5 crossings
    %% --------------------------------------------------------

    crossingX = [];

    for i = 1:NX-1

        f1 = phiLine(i)   - phiInterface;
        f2 = phiLine(i+1) - phiInterface;

        if f1*f2 <= 0

            if abs(phiLine(i+1)-phiLine(i)) > 1.0e-14

                fraction = ...
                    (phiInterface-phiLine(i)) ...
                    / ...
                    (phiLine(i+1)-phiLine(i));

                xc = ...
                    x(i) ...
                    + ...
                    fraction*(x(i+1)-x(i));

                crossingX(end+1) = xc;

            end
        end
    end

    numberOfCrossings(n) = length(crossingX);

    %% --------------------------------------------------------
    % Outermost left/right intersections
    %% --------------------------------------------------------

    xLeft = min(crossingX);

    xRight = max(crossingX);

    rNumerical(n) = ...
        0.5*(xRight-xLeft);

    rOverD(n) = ...
        rNumerical(n)/D;

    %% --------------------------------------------------------
    % Reference power law
    %
    % r/D = 1.35 sqrt(T)
    %% --------------------------------------------------------

    rTheoryOverD(n) = ...
        1.35*sqrt(T(n));

    relativeDifference(n) = ...
        abs( ...
            rOverD(n)-rTheoryOverD(n) ...
        ) ...
        / ...
        rTheoryOverD(n);

end

%% ============================================================
% Display numerical results
%% ============================================================

fprintf('\n');
fprintf('=============================================\n');
fprintf('Droplet impact spreading radius\n');
fprintf('Probe position: y = h + W = %.4f\n',yProbe);
fprintf('=============================================\n');

fprintf( ...
    '%8s %14s %14s %14s %14s %10s\n', ...
    'T', ...
    'r_num', ...
    'r_num/D', ...
    '1.35sqrt(T)', ...
    'rel.diff.', ...
    'Ncross' ...
);

for n = 1:length(T)

    fprintf( ...
        '%8.3f %14.6f %14.6f %14.6f %14.6e %10d\n', ...
        T(n), ...
        rNumerical(n), ...
        rOverD(n), ...
        rTheoryOverD(n), ...
        relativeDifference(n), ...
        numberOfCrossings(n) ...
    );

end

%% ============================================================
% Plot
%% ============================================================

figure;

loglog( ...
    T, ...
    rTheoryOverD, ...
    '-', ...
    'LineWidth', ...
    1.5 ...
);

hold on;

loglog( ...
    T, ...
    rOverD, ...
    'o', ...
    'MarkerSize', ...
    7, ...
    'LineWidth', ...
    1.2 ...
);

grid on;

xlabel('T = V_0 t / D');

ylabel('r / D');

legend( ...
    '1.35\sqrt{T}', ...
    'Present CHE', ...
    'Location', ...
    'northwest' ...
);

title('Droplet impact spreading radius');