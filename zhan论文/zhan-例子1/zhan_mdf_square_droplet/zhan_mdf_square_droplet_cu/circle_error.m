clear;
clc;
close all;

% ============================================================
% Final square-droplet interface error
%
% Numerical interface:
%       phi = 0
%
% Theoretical equilibrium:
%       R = 0.5
% ============================================================


% ------------------------------------------------------------
% 1. Load final phase field
% ------------------------------------------------------------

filename = 'phi_2000000.dat';

phi = load(filename);

[Ny, Nx] = size(phi);


% ------------------------------------------------------------
% 2. Physical coordinates
%
% Domain:
%       [-1, 1] x [-1, 1]
%
% Cell-centered coordinates
% ------------------------------------------------------------

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


% ------------------------------------------------------------
% 3. Extract phi = 0 contour
%
% contourc returns packed contour data.
% ------------------------------------------------------------

C = contourc( ...
    x, ...
    y, ...
    phi, ...
    [0 0]);


% ------------------------------------------------------------
% 4. Read every contour segment
%
% Usually there is only one closed droplet interface.
% To be safe, choose the segment containing the most points.
% ------------------------------------------------------------

k = 1;

bestN = 0;

xc = [];
yc = [];


while k < size(C, 2)

    level = C(1, k);

    nPoints = round(C(2, k));

    points = C(:, k+1 : k+nPoints);


    if nPoints > bestN

        bestN = nPoints;

        xc = points(1, :);

        yc = points(2, :);

    end


    k = k + nPoints + 1;

end


fprintf('Number of interface points = %d\n', bestN);


% ------------------------------------------------------------
% 5. Convert contour to polar coordinates
%
% Theoretical center = (0,0)
% ------------------------------------------------------------

theta = atan2(yc, xc);

theta = mod(theta, 2*pi);


r = sqrt( ...
    xc.^2 ...
    + ...
    yc.^2);


% ------------------------------------------------------------
% 6. Sort according to theta
% ------------------------------------------------------------

[thetaSort, order] = sort(theta);

rSort = r(order);


% Remove repeated theta values
[thetaUnique, ia] = unique( ...
    thetaSort, ...
    'stable');

rUnique = rSort(ia);

thetaUnique = thetaUnique(:);
rUnique = rUnique(:);

% ------------------------------------------------------------
% 7. Periodic extension
%
% Needed for interpolation near theta = 0 / 2*pi
% ------------------------------------------------------------

thetaExt = [ ...
    thetaUnique(end) - 2*pi; ...
    thetaUnique; ...
    thetaUnique(1) + 2*pi];


rExt = [ ...
    rUnique(end); ...
    rUnique; ...
    rUnique(1)];


% ------------------------------------------------------------
% 8. Uniform angular grid
%
% 720 directions
% ------------------------------------------------------------

thetaUniform = linspace( ...
    0, ...
    2*pi, ...
    721);

thetaUniform(end) = [];
thetaUniform = thetaUniform(:);

rUniform = interp1( ...
    thetaExt, ...
    rExt, ...
    thetaUniform, ...
    'linear');


% ------------------------------------------------------------
% 9. Compare with theoretical radius
% ------------------------------------------------------------

R = 0.5;

errorR = ...
    rUniform ...
    - ...
    R;


absError = abs(errorR);


% Maximum absolute radius error
maxError = max(absError);


% Mean absolute error
MAE = mean(absError);


% Root-mean-square error
RMSE = sqrt( ...
    mean( ...
        errorR.^2 ...
    ));


% Relative maximum error
maxRelativeError = ...
    maxError ...
    / ...
    R;


% Relative RMS error
relativeRMSE = ...
    RMSE ...
    / ...
    R;


% Mean numerical radius
meanRadius = mean(rUniform);


% Minimum / maximum radius
minRadius = min(rUniform);

maxRadius = max(rUniform);


% ------------------------------------------------------------
% 10. Area error
% ------------------------------------------------------------

areaNumerical = polyarea( ...
    xc, ...
    yc);


areaTheory = ...
    pi ...
    * ...
    R^2;


relativeAreaError = ...
    abs( ...
        areaNumerical ...
        - ...
        areaTheory ...
    ) ...
    / ...
    areaTheory;


% Equivalent radius calculated from area
equivalentRadius = sqrt( ...
    areaNumerical ...
    / ...
    pi);


% ------------------------------------------------------------
% 11. Print results
% ------------------------------------------------------------

fprintf('\n');
fprintf('=============================================\n');
fprintf('Final equilibrium circle error\n');
fprintf('=============================================\n');

fprintf('Theoretical radius R       = %.12e\n', R);

fprintf('Mean numerical radius      = %.12e\n', meanRadius);

fprintf('Minimum numerical radius   = %.12e\n', minRadius);

fprintf('Maximum numerical radius   = %.12e\n', maxRadius);

fprintf('\n');

fprintf('Max |r-R|                  = %.12e\n', maxError);

fprintf('MAE                        = %.12e\n', MAE);

fprintf('RMSE                       = %.12e\n', RMSE);

fprintf('Max relative error         = %.12e\n', maxRelativeError);

fprintf('Relative RMSE              = %.12e\n', relativeRMSE);

fprintf('\n');

fprintf('Numerical area             = %.12e\n', areaNumerical);

fprintf('Theoretical area           = %.12e\n', areaTheory);

fprintf('Relative area error        = %.12e\n', relativeAreaError);

fprintf('Equivalent radius          = %.12e\n', equivalentRadius);

fprintf('=============================================\n');


% ------------------------------------------------------------
% 12. Plot numerical radius against theory
% ------------------------------------------------------------

figure;

plot( ...
    thetaUniform, ...
    rUniform, ...
    'LineWidth', ...
    1.6);

hold on;


yline( ...
    R, ...
    'k--', ...
    'LineWidth', ...
    1.5);


xlabel('\theta');

ylabel('r(\theta)');

title('Final interface radius');

legend( ...
    'Numerical', ...
    'R = 0.5', ...
    'Location', ...
    'best');

grid on;

box on;


% ------------------------------------------------------------
% 13. Plot radius error
% ------------------------------------------------------------

figure;

plot( ...
    thetaUniform, ...
    errorR, ...
    'LineWidth', ...
    1.5);


yline( ...
    0, ...
    'k--');


xlabel('\theta');

ylabel('r(\theta) - 0.5');

title('Final interface radius error');

grid on;

box on;