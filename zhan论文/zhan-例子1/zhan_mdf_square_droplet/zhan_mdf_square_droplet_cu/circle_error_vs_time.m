clear;
clc;
close all;

% ============================================================
% Square droplet:
% error of phi = 0 interface against theoretical circle R = 0.5
%
% Files:
% phi_00000.dat
% phi_100000.dat
% ...
% phi_2000000.dat
% ============================================================


% ------------------------------------------------------------
% 1. Parameters
% ------------------------------------------------------------

steps = 0:100000:2000000;

dt = 1.0e-4;

R = 0.5;

nTheta = 720;

nCase = length(steps);


% ------------------------------------------------------------
% 2. Arrays for results
% ------------------------------------------------------------

physicalTime = zeros(nCase,1);

meanRadius   = zeros(nCase,1);
minRadius    = zeros(nCase,1);
maxRadius    = zeros(nCase,1);

maxError     = zeros(nCase,1);
MAE          = zeros(nCase,1);
RMSE         = zeros(nCase,1);

relativeRMSE = zeros(nCase,1);

areaNum      = zeros(nCase,1);
areaError    = zeros(nCase,1);

Req          = zeros(nCase,1);


% ------------------------------------------------------------
% 3. Loop over all files
% ------------------------------------------------------------

for n = 1:nCase

    step = steps(n);

    physicalTime(n) = step * dt;

    filename = sprintf('phi_%05d.dat', step);


    if ~isfile(filename)

        warning('Cannot find %s', filename);

        meanRadius(n)   = NaN;
        minRadius(n)    = NaN;
        maxRadius(n)    = NaN;
        maxError(n)     = NaN;
        MAE(n)          = NaN;
        RMSE(n)         = NaN;
        relativeRMSE(n) = NaN;
        areaNum(n)      = NaN;
        areaError(n)    = NaN;
        Req(n)          = NaN;

        continue;

    end


    % --------------------------------------------------------
    % Load phi
    % --------------------------------------------------------

    phi = load(filename);

    [Ny, Nx] = size(phi);


    dx = 2.0 / Nx;
    dy = 2.0 / Ny;


    x = linspace( ...
        -1.0 + dx/2, ...
         1.0 - dx/2, ...
         Nx);


    y = linspace( ...
        -1.0 + dy/2, ...
         1.0 - dy/2, ...
         Ny);


    % --------------------------------------------------------
    % Extract phi = 0 contour
    % --------------------------------------------------------

    C = contourc( ...
        x, ...
        y, ...
        phi, ...
        [0 0]);


    % --------------------------------------------------------
    % Find largest contour segment
    % --------------------------------------------------------

    k = 1;

    bestN = 0;

    xc = [];
    yc = [];


    while k < size(C,2)

        nPoints = round(C(2,k));

        points = C(:, k+1:k+nPoints);


        if nPoints > bestN

            bestN = nPoints;

            xc = points(1,:);
            yc = points(2,:);

        end


        k = k + nPoints + 1;

    end


    % --------------------------------------------------------
    % Polar coordinates
    % --------------------------------------------------------

    theta = atan2(yc, xc);

    theta = mod(theta, 2*pi);


    r = sqrt( ...
        xc.^2 ...
        + ...
        yc.^2);


    % --------------------------------------------------------
    % Sort by theta
    % --------------------------------------------------------

    [thetaSort, order] = sort(theta);

    rSort = r(order);


    [thetaUnique, ia] = unique( ...
        thetaSort, ...
        'stable');


    rUnique = rSort(ia);


    % Force column vectors
    thetaUnique = thetaUnique(:);

    rUnique = rUnique(:);


    % --------------------------------------------------------
    % Periodic extension
    % --------------------------------------------------------

    thetaExt = [ ...
        thetaUnique(end) - 2*pi; ...
        thetaUnique; ...
        thetaUnique(1) + 2*pi];


    rExt = [ ...
        rUnique(end); ...
        rUnique; ...
        rUnique(1)];


    % --------------------------------------------------------
    % Uniform theta
    % --------------------------------------------------------

    thetaUniform = linspace( ...
        0, ...
        2*pi, ...
        nTheta + 1);


    thetaUniform(end) = [];

    thetaUniform = thetaUniform(:);


    rUniform = interp1( ...
        thetaExt, ...
        rExt, ...
        thetaUniform, ...
        'linear');


    % --------------------------------------------------------
    % Radius errors
    % --------------------------------------------------------

    e = rUniform - R;


    meanRadius(n) = mean(rUniform);

    minRadius(n) = min(rUniform);

    maxRadius(n) = max(rUniform);


    maxError(n) = max(abs(e));

    MAE(n) = mean(abs(e));

    RMSE(n) = sqrt(mean(e.^2));


    relativeRMSE(n) = RMSE(n) / R;


    % --------------------------------------------------------
    % Area
    % --------------------------------------------------------

    areaNum(n) = polyarea(xc, yc);


    areaTheory = pi * R^2;


    areaError(n) = ...
        abs(areaNum(n) - areaTheory) ...
        / ...
        areaTheory;


    Req(n) = sqrt(areaNum(n) / pi);


    % --------------------------------------------------------
    % Print each step
    % --------------------------------------------------------

    fprintf( ...
        't = %6.1f   meanR = %.8f   RMSErel = %.6e   AreaErr = %.6e\n', ...
        physicalTime(n), ...
        meanRadius(n), ...
        relativeRMSE(n), ...
        areaError(n));

end


% ============================================================
% 4. Result table
% ============================================================

Result = table( ...
    steps(:), ...
    physicalTime, ...
    meanRadius, ...
    minRadius, ...
    maxRadius, ...
    maxError, ...
    RMSE, ...
    relativeRMSE, ...
    areaError, ...
    Req, ...
    'VariableNames', ...
    { ...
    'Step', ...
    'Time', ...
    'MeanRadius', ...
    'MinRadius', ...
    'MaxRadius', ...
    'MaxError', ...
    'RMSE', ...
    'RelativeRMSE', ...
    'RelativeAreaError', ...
    'EquivalentRadius' ...
    } ...
    );


disp(Result);


% ============================================================
% 5. Plot relative RMS error versus time
% ============================================================

figure;

plot( ...
    physicalTime, ...
    100 * relativeRMSE, ...
    'o-', ...
    'LineWidth', 1.5);

xlabel('t');

ylabel('Relative RMSE (%)');

title('Radius error versus time');

grid on;
box on;


% ============================================================
% 6. Plot equivalent radius
% ============================================================

figure;

plot( ...
    physicalTime, ...
    Req, ...
    'o-', ...
    'LineWidth', 1.5);

hold on;


yline( ...
    R, ...
    'k--', ...
    'LineWidth', 1.5);


xlabel('t');

ylabel('Equivalent radius');

title('Equivalent radius versus time');

legend( ...
    'Numerical', ...
    'R = 0.5', ...
    'Location', ...
    'best');

grid on;
box on;


% ============================================================
% 7. Plot area error
% ============================================================

figure;

plot( ...
    physicalTime, ...
    100 * areaError, ...
    'o-', ...
    'LineWidth', 1.5);

xlabel('t');

ylabel('Relative area error (%)');

title('Area error versus time');

grid on;
box on;