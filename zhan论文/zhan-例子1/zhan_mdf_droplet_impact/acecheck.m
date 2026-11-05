clear;
clc;
close all;

%% ============================================================
%  Zhan 2023 droplet impact
%  ACE, rhoA:rhoB = 1000:1
%
%  Spreading radius post-processing
%
%  Definition:
%
%  At y = h + W,
%  find all phi = 0.5 intersections.
%
%  Use the outermost left/right intersections:
%
%        r = (xR - xL)/2
%
%  Compare with:
%
%        r/D = 1.35 sqrt(T)
%
%  IMPORTANT:
%  This is our reproducible post-processing definition.
%% ============================================================


%% ============================================================
% 1. ACE files
%% ============================================================

files = {
    'droplet_phi_ACE_iso_step_00400.dat'
    'droplet_phi_ACE_iso_step_00800.dat'
    'droplet_phi_ACE_iso_step_01600.dat'
    'droplet_phi_ACE_iso_step_03200.dat'
    'droplet_phi_ACE_iso_step_06400.dat'
};


%% ============================================================
% 2. Dimensionless times
%% ============================================================

T = [
    0.1
    0.2
    0.4
    0.8
    1.6
];


%% ============================================================
% 3. Parameters
%% ============================================================

NX = 1500;
NY = 500;

dx = 1.0;

D = 200.0;

h = 50.0;

W = 5.0;

phiInterface = 0.5;


% Probe position
%
% y = h + W = 55
%
yProbe = h + W;


%% ============================================================
% 4. Coordinates
%% ============================================================

x = -750 + (0.5 : 1 : 1499.5);

y = 0.5 : 1 : 499.5;


%% ============================================================
% 5. Output arrays
%% ============================================================

rNumerical = zeros(length(T),1);

rOverD = zeros(length(T),1);

rTheoryOverD = zeros(length(T),1);

relativeDifference = zeros(length(T),1);

numberOfCrossings = zeros(length(T),1);

xLeftAll = zeros(length(T),1);

xRightAll = zeros(length(T),1);


%% ============================================================
% 6. Loop over ACE files
%% ============================================================

for n = 1:length(T)

    fprintf('\nReading: %s\n',files{n});

    phi = readmatrix(files{n});


    %% --------------------------------------------------------
    % Check matrix size
    %% --------------------------------------------------------

    if size(phi,1) ~= NY || size(phi,2) ~= NX

        error( ...
            'Unexpected matrix size in %s. Got %d x %d.', ...
            files{n}, ...
            size(phi,1), ...
            size(phi,2) ...
        );

    end


    %% --------------------------------------------------------
    % Interpolate phi onto y = yProbe
    %
    % Grid centers are:
    %
    % y = 0.5, 1.5, 2.5, ...
    %
    % Therefore y = 55 lies between two grid lines.
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

        f1 = ...
            phiLine(i) ...
            - ...
            phiInterface;


        f2 = ...
            phiLine(i+1) ...
            - ...
            phiInterface;


        % True crossing
        if f1*f2 <= 0

            denominator = ...
                phiLine(i+1) ...
                - ...
                phiLine(i);


            if abs(denominator) > 1.0e-14

                fraction = ...
                    ( ...
                        phiInterface ...
                        - ...
                        phiLine(i) ...
                    ) ...
                    / ...
                    denominator;


                xc = ...
                    x(i) ...
                    + ...
                    fraction ...
                    * ...
                    (x(i+1)-x(i));


                crossingX(end+1) = xc;

            end

        end

    end


    %% --------------------------------------------------------
    % Number of crossings
    %% --------------------------------------------------------

    numberOfCrossings(n) = length(crossingX);


    if isempty(crossingX)

        error( ...
            'No phi = 0.5 crossing found at T = %.3f.', ...
            T(n) ...
        );

    end


    %% --------------------------------------------------------
    % Outermost left/right intersections
    %% --------------------------------------------------------

    xLeft = min(crossingX);

    xRight = max(crossingX);


    xLeftAll(n) = xLeft;

    xRightAll(n) = xRight;


    %% --------------------------------------------------------
    % Numerical spreading radius
    %% --------------------------------------------------------

    rNumerical(n) = ...
        0.5 ...
        * ...
        (xRight-xLeft);


    rOverD(n) = ...
        rNumerical(n) ...
        / ...
        D;


    %% --------------------------------------------------------
    % Reference power law
    %
    % r/D = 1.35 sqrt(T)
    %% --------------------------------------------------------

    rTheoryOverD(n) = ...
        1.35 ...
        * ...
        sqrt(T(n));


    relativeDifference(n) = ...
        abs( ...
            rOverD(n) ...
            - ...
            rTheoryOverD(n) ...
        ) ...
        / ...
        rTheoryOverD(n);

end


%% ============================================================
% 7. Display numerical results
%% ============================================================

fprintf('\n');
fprintf('===============================================================\n');
fprintf('Zhan 2023 droplet impact - ACE\n');
fprintf('rhoA : rhoB = 1000 : 1\n');
fprintf('Probe position: y = h + W = %.4f\n',yProbe);
fprintf('===============================================================\n');

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
% 8. Also display outer intersections
%% ============================================================

fprintf('\n');
fprintf('Outer phi = 0.5 intersections\n');
fprintf('===============================================================\n');

fprintf( ...
    '%8s %16s %16s\n', ...
    'T', ...
    'xLeft', ...
    'xRight' ...
);


for n = 1:length(T)

    fprintf( ...
        '%8.3f %16.6f %16.6f\n', ...
        T(n), ...
        xLeftAll(n), ...
        xRightAll(n) ...
    );

end


%% ============================================================
% 9. Plot spreading radius
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


xlabel( ...
    'T = V_0 t / D' ...
);


ylabel( ...
    'r / D' ...
);


legend( ...
    '1.35\sqrt{T}', ...
    'Present ACE, \rho_A/\rho_B = 1000', ...
    'Location', ...
    'northwest' ...
);


title( ...
    'Droplet impact spreading radius - ACE' ...
);


%% ============================================================
% 10. Plot relative difference
%% ============================================================

figure;

semilogx( ...
    T, ...
    100*relativeDifference, ...
    'o-', ...
    'LineWidth', ...
    1.2, ...
    'MarkerSize', ...
    7 ...
);


grid on;


xlabel( ...
    'T = V_0 t / D' ...
);


ylabel( ...
    'Relative difference (%)' ...
);


title( ...
    'Difference from r/D = 1.35\sqrt{T}' ...
);


%% ============================================================
% 11. Effective coefficient
%
% r/D = C sqrt(T)
%% ============================================================

Cnumerical = ...
    rOverD ...
    ./ ...
    sqrt(T);


fprintf('\n');
fprintf('===============================================================\n');
fprintf('Effective spreading coefficient\n');
fprintf('r/D = C sqrt(T)\n');
fprintf('===============================================================\n');

fprintf( ...
    '%8s %14s\n', ...
    'T', ...
    'C_num' ...
);


for n = 1:length(T)

    fprintf( ...
        '%8.3f %14.6f\n', ...
        T(n), ...
        Cnumerical(n) ...
    );

end


fprintf('\nReference coefficient C = 1.35\n');