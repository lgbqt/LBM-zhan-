clear;
clc;
close all;

%% ============================================================
% Zhan 2023 droplet impact
% ACE, rhoA/rhoB = 1000
% Re = 500
%
% Spreading-radius post-processing
%
% Definition used here:
%
%   yProbe = h + W = 55
%
% At this height find phi = 0.5 crossings and use
% outermost left/right intersections:
%
%       r = (xRight-xLeft)/2
%% ============================================================

files = {
    'droplet_phi_ACE_iso_step_00400.dat'
    'droplet_phi_ACE_iso_step_00800.dat'
    'droplet_phi_ACE_iso_step_01600.dat'
    'droplet_phi_ACE_iso_step_03200.dat'
    'droplet_phi_ACE_iso_step_06400.dat'
};

T = [
    0.1
    0.2
    0.4
    0.8
    1.6
];

%% Parameters

NX = 1500;
NY = 500;

D = 200.0;

h = 50.0;
W = 5.0;

phiInterface = 0.5;

yProbe = h + W;

x = -750 + (0.5:1:1499.5);
y = 0.5:1:499.5;

%% Results

rNumerical = zeros(length(T),1);
rOverD = zeros(length(T),1);
rTheory = 1.35*sqrt(T);

relativeDifference = zeros(length(T),1);
numberOfCrossings = zeros(length(T),1);

xLeftAll = zeros(length(T),1);
xRightAll = zeros(length(T),1);

%% ============================================================
% Loop
%% ============================================================

for n = 1:length(T)

    phi = readmatrix(files{n});

    if size(phi,1) ~= NY || size(phi,2) ~= NX
        error('Unexpected data size in %s',files{n});
    end

    %% Interpolate onto y = 55

    upperIndex = find(y > yProbe,1);
    lowerIndex = upperIndex - 1;

    y1 = y(lowerIndex);
    y2 = y(upperIndex);

    alpha = ...
        (yProbe-y1)/(y2-y1);

    phiLine = ...
        (1-alpha)*phi(lowerIndex,:) ...
        + alpha*phi(upperIndex,:);

    %% Find phi = 0.5 crossings

    crossingX = [];

    for i = 1:NX-1

        f1 = phiLine(i)   - phiInterface;
        f2 = phiLine(i+1) - phiInterface;

        if f1*f2 <= 0

            denominator = ...
                phiLine(i+1)-phiLine(i);

            if abs(denominator) > 1.0e-14

                fraction = ...
                    (phiInterface-phiLine(i)) ...
                    / denominator;

                xc = ...
                    x(i) ...
                    + fraction*(x(i+1)-x(i));

                crossingX(end+1) = xc;
            end
        end
    end

    numberOfCrossings(n) = length(crossingX);

    if isempty(crossingX)
        error('No phi=0.5 crossing at T = %.3f',T(n));
    end

    %% Outermost intersections

    xLeft = min(crossingX);
    xRight = max(crossingX);

    xLeftAll(n) = xLeft;
    xRightAll(n) = xRight;

    rNumerical(n) = ...
        0.5*(xRight-xLeft);

    rOverD(n) = ...
        rNumerical(n)/D;

    relativeDifference(n) = ...
        abs(rOverD(n)-rTheory(n)) ...
        / rTheory(n);

end

%% ============================================================
% Output
%% ============================================================

Cnum = ...
    rOverD ./ sqrt(T);

fprintf('\n');
fprintf('===============================================================\n');
fprintf('Zhan 2023 droplet impact - ACE\n');
fprintf('rhoA/rhoB = 1000, Re = 500\n');
fprintf('Probe position y = %.4f\n',yProbe);
fprintf('===============================================================\n');

fprintf( ...
    '%8s %14s %14s %14s %14s %10s %12s\n', ...
    'T', ...
    'r_num', ...
    'r/D', ...
    '1.35sqrt(T)', ...
    'rel.diff.', ...
    'Ncross', ...
    'C_num');

for n = 1:length(T)

    fprintf( ...
        '%8.3f %14.6f %14.6f %14.6f %14.6e %10d %12.6f\n', ...
        T(n), ...
        rNumerical(n), ...
        rOverD(n), ...
        rTheory(n), ...
        relativeDifference(n), ...
        numberOfCrossings(n), ...
        Cnum(n));

end

%% ============================================================
% Plot
%% ============================================================

figure;

loglog( ...
    T, ...
    rTheory, ...
    '-', ...
    'LineWidth',1.5);

hold on;

loglog( ...
    T, ...
    rOverD, ...
    'o', ...
    'MarkerSize',7, ...
    'LineWidth',1.2);

grid on;

xlabel('T = V_0 t / D');
ylabel('r / D');

legend( ...
    '1.35\sqrt{T}', ...
    'Present ACE, Re = 500', ...
    'Location','northwest');

title('Droplet impact spreading radius, Re = 500');