clear;
clc;
close all;

% ============================================================
% Zhan 2023 layered Poiseuille flow
% Error-distribution analysis
%
% numerical file columns:
% 1 y
% 2 phi
% 3 rho
% 4 u_num
% 5 u_exact
% 6 du_num
% 7 du_exact
% ============================================================


% ------------------------------------------------------------
% 1. Load data
% ------------------------------------------------------------

filename = 'poiseuille_numerical_100to1_B.dat';

data = load(filename);

y       = data(:,1);
phi     = data(:,2);
rho     = data(:,3);

uNum    = data(:,4);
uExact  = data(:,5);

duNum   = data(:,6);
duExact = data(:,7);


% ------------------------------------------------------------
% 2. Pointwise errors
% ------------------------------------------------------------

errorU = uNum - uExact;

absErrorU = abs(errorU);

errorDu = duNum - duExact;

absErrorDu = abs(errorDu);


% ------------------------------------------------------------
% 3. Paper Eq.(52) errors
% ------------------------------------------------------------

ErrU = sum(absErrorU) / sum(abs(uExact));

ErrDu = sum(absErrorDu) / sum(abs(duExact));


fprintf('=============================================\n');
fprintf('Poiseuille error analysis\n');
fprintf('=============================================\n');

fprintf('Err_u       = %.12e\n', ErrU);
fprintf('Err_du/dy   = %.12e\n', ErrDu);

fprintf('\n');


% ------------------------------------------------------------
% 4. Find largest pointwise velocity error
% ------------------------------------------------------------

[maxAbsU, indexU] = max(absErrorU);

fprintf('Maximum pointwise |u error|\n');
fprintf('value       = %.12e\n', maxAbsU);
fprintf('y           = %.12e\n', y(indexU));
fprintf('u_num       = %.12e\n', uNum(indexU));
fprintf('u_exact     = %.12e\n', uExact(indexU));

fprintf('\n');


% ------------------------------------------------------------
% 5. Find largest pointwise gradient error
% ------------------------------------------------------------

[maxAbsDu, indexDu] = max(absErrorDu);

fprintf('Maximum pointwise |du/dy error|\n');
fprintf('value       = %.12e\n', maxAbsDu);
fprintf('y           = %.12e\n', y(indexDu));
fprintf('du_num      = %.12e\n', duNum(indexDu));
fprintf('du_exact    = %.12e\n', duExact(indexDu));

fprintf('=============================================\n');


% ============================================================
% 6. Velocity comparison
% ============================================================

figure;

plot(uExact, y, 'k--', 'LineWidth', 2);
hold on;

plot(uNum, y, 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);

xlabel('u');
ylabel('y');

title('Layered Poiseuille velocity');

legend('Analytical', 'Numerical', 'Location', 'best');

grid on;
box on;


% ============================================================
% 7. Velocity error versus y
% ============================================================

figure;

plot(errorU, y, 'LineWidth', 1.6);

xline(0, 'k--');

xlabel('u_{num} - u_{exact}');
ylabel('y');

title('Velocity error distribution');

grid on;
box on;


% ============================================================
% 8. Absolute velocity error versus y
% ============================================================

figure;

plot(absErrorU, y, 'LineWidth', 1.6);

xlabel('|u_{num} - u_{exact}|');
ylabel('y');

title('Absolute velocity error distribution');

grid on;
box on;


% ============================================================
% 9. du/dy comparison
% ============================================================

figure;

plot(duExact, y, 'k--', 'LineWidth', 2);
hold on;

plot(duNum, y, 'o-', 'LineWidth', 1.2, 'MarkerSize', 3);

xlabel('\partial u / \partial y');
ylabel('y');

title('Velocity-gradient comparison');

legend('Analytical', 'Numerical', 'Location', 'best');

grid on;
box on;


% ============================================================
% 10. du/dy error versus y
% ============================================================

figure;

plot(errorDu, y, 'LineWidth', 1.6);

xline(0, 'k--');

xlabel('(du/dy)_{num} - (du/dy)_{exact}');
ylabel('y');

title('Velocity-gradient error distribution');

grid on;
box on;


% ============================================================
% 11. Absolute du/dy error
% ============================================================

figure;

plot(absErrorDu, y, 'LineWidth', 1.6);

xlabel('|(du/dy)_{num} - (du/dy)_{exact}|');
ylabel('y');

title('Absolute velocity-gradient error');

grid on;
box on;


% ============================================================
% 12. Also plot phi and rho
%
% This helps us identify whether the error peak occurs
% near the diffuse interface y = 0.
% ============================================================

figure;

plot(phi, y, 'LineWidth', 1.6);

xlabel('\phi');
ylabel('y');

title('Final phase-field profile');

grid on;
box on;


figure;

plot(rho, y, 'LineWidth', 1.6);

xlabel('\rho');
ylabel('y');

title('Final density profile');

grid on;
box on;