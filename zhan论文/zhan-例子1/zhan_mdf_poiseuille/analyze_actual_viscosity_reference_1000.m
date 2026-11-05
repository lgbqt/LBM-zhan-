
clc;
close all;

% ============================================================
% Zhan 2023 layered Poiseuille flow
%
% Diagnostic:
% Compare
%
% 1. LBM numerical velocity
% 2. Paper analytical velocity based on Eq.(49) tanh profile
% 3. New reference velocity based on FINAL numerical viscosity
%
% Current case:
% rhoA : rhoB = 1000 : 1
% ============================================================


% ============================================================
% 1. Load final LBM data
% ============================================================

filename = 'poiseuille_numerical_1000to1_A.dat';

data = load(filename);

y        = data(:,1);
phiNum   = data(:,2);
rhoNum   = data(:,3);

uLBM     = data(:,4);
uPaper   = data(:,5);

duLBM    = data(:,6);
duPaper  = data(:,7);


% ============================================================
% 2. Parameters
% ============================================================

phiA = 1.0;
phiB = -1.0;

rhoA = 1000.0;
rhoB = 1.0;

nuA = 0.1;
nuB = 0.1;

W = 0.04;

h = 0.5;

Uc = 1.0e-4;


% This is the Gx already obtained from Eq.(50)-(51)
% using the paper tanh interface.

Gx = 3.708100368381e-02;

% ============================================================
% 3. Original Eq.(49) phase profile
% ============================================================

phiExact = ...
    0.5 * (phiA + phiB) ...
    + 0.5 * (phiA - phiB) .* tanh(-2.0 .* y ./ W);


% ============================================================
% 4. Original analytical density / dynamic viscosity
% ============================================================

rhoExact = ...
    (rhoA - rhoB) ...
    / (phiA - phiB) ...
    .* (phiExact - phiB) ...
    + rhoB;


% Here nuA = nuB = 0.1
%
% Dynamic viscosity:
%
% mu = rho * nu
%
% More generally the paper interpolates dynamic viscosity.

muExact = rhoExact .* nuA;


% ============================================================
% 5. FINAL numerical dynamic viscosity
% ============================================================

muNum = rhoNum .* nuA;


% ============================================================
% 6. Difference between final interface and Eq.(49)
% ============================================================

deltaPhi = phiNum - phiExact;

deltaRho = rhoNum - rhoExact;

deltaMu = muNum - muExact;


fprintf('=============================================\n');
fprintf('Final-interface diagnostic\n');
fprintf('=============================================\n');

fprintf('max |phi_num - phi_exact| = %.12e\n', ...
    max(abs(deltaPhi)));

fprintf('max |rho_num - rho_exact| = %.12e\n', ...
    max(abs(deltaRho)));

fprintf('max |mu_num - mu_exact|   = %.12e\n', ...
    max(abs(deltaMu)));

fprintf('\n');


% ============================================================
% 7. Build a fine y grid
%
% We interpolate the FINAL numerical viscosity onto this grid.
%
% Fluid nodes are at:
%
% -0.495, -0.485, ..., 0.495
%
% Walls are at:
%
% -0.5 and 0.5
%
% At the two very small boundary gaps we simply extend the
% nearest viscosity value.
% ============================================================

NFINE = 20001;

yFine = linspace(-h, h, NFINE)';


% Add wall points for interpolation.
%
% mu(-0.5) = first fluid-cell value
% mu(+0.5) = last fluid-cell value

yExtended = [-h; y; h];

muExtended = [muNum(1); muNum; muNum(end)];


muFine = interp1( ...
    yExtended, ...
    muExtended, ...
    yFine, ...
    'pchip' ...
);


% ============================================================
% 8. New reference based on FINAL numerical viscosity
%
% Steady equation:
%
% d/dy( mu du/dy ) + Gx = 0
%
% Integrate once:
%
% mu du/dy = -Gx*y + C1
%
% Boundary:
%
% u(-h)=0
% u(+h)=0
%
% Define
%
% A_num = integral y/mu dy
% B_num = integral 1/mu dy
%
% Then
%
% C1 = Gx * A_num / B_num
% ============================================================

integrandA = yFine ./ muFine;

integrandB = 1.0 ./ muFine;


A_num = trapz(yFine, integrandA);

B_num = trapz(yFine, integrandB);


C1_num = Gx * A_num / B_num;


fprintf('Reference using FINAL numerical viscosity\n');

fprintf('A_num                  = %.12e\n', A_num);
fprintf('B_num                  = %.12e\n', B_num);
fprintf('C1_num                 = %.12e\n', C1_num);

fprintf('\n');


% ============================================================
% 9. Compute new reference du/dy
% ============================================================

duActualFine = ...
    (-Gx .* yFine + C1_num) ...
    ./ muFine;


% ============================================================
% 10. Integrate du/dy to obtain u
%
% u(-h) = 0
% ============================================================

uActualFine = cumtrapz(yFine, duActualFine);


% Because of numerical integration, remove any tiny linear
% end-point residual so that:
%
% u(-h)=u(+h)=0
%
% exactly on the fine grid.

endResidual = uActualFine(end);

linearCorrection = ...
    endResidual ...
    .* (yFine + h) ...
    ./ (2.0 * h);

uActualFine = ...
    uActualFine ...
    - linearCorrection;


% ============================================================
% 11. Interpolate new reference back to LBM cell centers
% ============================================================

uActual = interp1( ...
    yFine, ...
    uActualFine, ...
    y, ...
    'pchip' ...
);


duActual = interp1( ...
    yFine, ...
    duActualFine, ...
    y, ...
    'pchip' ...
);


% ============================================================
% 12. Compute interface velocity of this new reference
%
% This is NOT forced to equal Uc anymore because Gx remains
% fixed at the original paper value.
% ============================================================

uActualAt0 = interp1( ...
    yFine, ...
    uActualFine, ...
    0.0, ...
    'pchip' ...
);


fprintf('u_actual_mu(0)          = %.12e\n', uActualAt0);
fprintf('original Uc             = %.12e\n', Uc);
fprintf('difference              = %.12e\n', uActualAt0 - Uc);

fprintf('\n');


% ============================================================
% 13. Eq.(52) error:
%
% LBM vs original paper analytical solution
% ============================================================

ErrPaper = ...
    sum(abs(uLBM - uPaper)) ...
    / sum(abs(uPaper));


ErrDuPaper = ...
    sum(abs(duLBM - duPaper)) ...
    / sum(abs(duPaper));


% ============================================================
% 14. Error:
%
% LBM vs reference using FINAL numerical viscosity
% ============================================================

ErrActual = ...
    sum(abs(uLBM - uActual)) ...
    / sum(abs(uActual));


ErrDuActual = ...
    sum(abs(duLBM - duActual)) ...
    / sum(abs(duActual));


fprintf('=============================================\n');
fprintf('Velocity errors\n');
fprintf('=============================================\n');

fprintf('LBM vs paper analytic\n');
fprintf('Err_u       = %.12e\n', ErrPaper);
fprintf('Err_du/dy   = %.12e\n', ErrDuPaper);

fprintf('\n');

fprintf('LBM vs FINAL-viscosity reference\n');
fprintf('Err_u       = %.12e\n', ErrActual);
fprintf('Err_du/dy   = %.12e\n', ErrDuActual);

fprintf('=============================================\n');


% ============================================================
% 15. Maximum errors and their positions
% ============================================================

errorPaper = uLBM - uPaper;

errorActual = uLBM - uActual;


[maxPaper, iPaper] = max(abs(errorPaper));

[maxActual, iActual] = max(abs(errorActual));


fprintf('\n');

fprintf('Maximum velocity error: original reference\n');
fprintf('|error|max = %.12e\n', maxPaper);
fprintf('y           = %.12e\n', y(iPaper));

fprintf('\n');

fprintf('Maximum velocity error: FINAL-mu reference\n');
fprintf('|error|max = %.12e\n', maxActual);
fprintf('y           = %.12e\n', y(iActual));

fprintf('\n');


% ============================================================
% 16. Plot phase field
% ============================================================

figure;

plot(phiExact, y, 'k--', 'LineWidth', 2);
hold on;

plot(phiNum, y, 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);

xlabel('\phi');
ylabel('y');

title('Phase field: Eq.(49) vs final LBM');

legend( ...
    'Eq.(49)', ...
    'Final LBM', ...
    'Location', ...
    'best' ...
);

grid on;
box on;


% ============================================================
% 17. Phase-field difference
% ============================================================

figure;

plot(deltaPhi, y, 'LineWidth', 1.6);

xline(0, 'k--');

xlabel('\phi_{final} - \phi_{Eq.(49)}');
ylabel('y');

title('Phase-field difference');

grid on;
box on;


% ============================================================
% 18. Dynamic-viscosity comparison
% ============================================================

figure;

plot(muExact, y, 'k--', 'LineWidth', 2);
hold on;

plot(muNum, y, 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);

xlabel('\mu');
ylabel('y');

title('Dynamic viscosity');

legend( ...
    'Eq.(49)-based', ...
    'Final LBM', ...
    'Location', ...
    'best' ...
);

grid on;
box on;


% ============================================================
% 19. Dynamic-viscosity difference
% ============================================================

figure;

plot(deltaMu, y, 'LineWidth', 1.6);

xline(0, 'k--');

xlabel('\mu_{final} - \mu_{Eq.(49)}');
ylabel('y');

title('Dynamic-viscosity difference');

grid on;
box on;


% ============================================================
% 20. Three velocity profiles
% ============================================================

figure;

plot(uPaper, y, 'k--', 'LineWidth', 2);
hold on;

plot(uActual, y, 'LineWidth', 1.8);

plot(uLBM, y, 'o-', 'LineWidth', 1.1, 'MarkerSize', 3);

xlabel('u');
ylabel('y');

title('Poiseuille velocity comparison');

legend( ...
    'Paper analytical', ...
    'Reference using final \mu', ...
    'LBM', ...
    'Location', ...
    'best' ...
);

grid on;
box on;


% ============================================================
% 21. Compare the two velocity errors
% ============================================================

figure;

plot(errorPaper, y, 'LineWidth', 1.6);
hold on;

plot(errorActual, y, 'LineWidth', 1.6);

xline(0, 'k--');

xlabel('Velocity error');
ylabel('y');

title('Velocity-error comparison');

legend( ...
    'LBM - paper analytical', ...
    'LBM - final-\mu reference', ...
    'Location', ...
    'best' ...
);

grid on;
box on;


% ============================================================
% 22. Gradient comparison
% ============================================================

figure;

plot(duPaper, y, 'k--', 'LineWidth', 2);
hold on;

plot(duActual, y, 'LineWidth', 1.8);

plot(duLBM, y, 'o-', 'LineWidth', 1.1, 'MarkerSize', 3);

xlabel('\partial u / \partial y');
ylabel('y');

title('Velocity-gradient comparison');

legend( ...
    'Paper analytical', ...
    'Reference using final \mu', ...
    'LBM local Eq.(40)', ...
    'Location', ...
    'best' ...
);

grid on;
box on;