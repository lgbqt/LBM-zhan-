clear;
clc;
close all;

% ============================================================
% Zhan 2023 square droplet
% GPU result: phi = 0 contour
%
% step = 0 : 100000 : 2000000
% dt   = 1e-4
%
% Therefore:
% t = 0, 10, 20, ..., 200
% ============================================================


% ------------------------------------------------------------
% 1. Time steps
% ------------------------------------------------------------

steps = 0:100000:2000000;

dt = 1.0e-4;


% ------------------------------------------------------------
% 2. Colors
% ------------------------------------------------------------

colors = lines(length(steps));


% ------------------------------------------------------------
% 3. Figure
% ------------------------------------------------------------

figure;
hold on;


% ------------------------------------------------------------
% 4. Read every phi file and draw phi = 0 contour
% ------------------------------------------------------------

for n = 1:length(steps)

    step = steps(n);

    filename = sprintf('phi_%05d.dat', step);

    if ~isfile(filename)

        warning('Cannot find %s', filename);

        continue;

    end


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


    t = step * dt;


    [~, h] = contour( ...
        x, ...
        y, ...
        phi, ...
        [0 0], ...
        'LineWidth', 1.5, ...
        'LineColor', colors(n,:) );


    h.DisplayName = sprintf('t = %.0f', t);

end

% ------------------------------------------------------------
% 5. Theoretical final circle
%
% Paper:
% final equilibrium radius R = 0.5
% ------------------------------------------------------------

theta = linspace(0, 2*pi, 500);

R = 0.5;

xCircle = R * cos(theta);
yCircle = R * sin(theta);


plot( ...
    xCircle, ...
    yCircle, ...
    'k--', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'R = 0.5');


% ------------------------------------------------------------
% 6. Figure settings
% ------------------------------------------------------------

axis equal;

xlim([-0.65 0.65]);
ylim([-0.65 0.65]);


xlabel('x');
ylabel('y');


title('Square droplet evolution (\phi = 0 contour)');


legend( ...
    'show', ...
    'Location', ...
    'eastoutside');


grid on;
box on;

hold off;