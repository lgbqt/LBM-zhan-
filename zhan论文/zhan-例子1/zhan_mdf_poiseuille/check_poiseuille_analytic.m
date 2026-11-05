clear;
clc;
close all;

% ============================================================
% Zhan 2023 layered Poiseuille flow
% Analytical solution check
% First case: rhoA : rhoB = 100 : 1
% ============================================================

% ------------------------------------------------------------
% 1. Parameters
% ------------------------------------------------------------

h = 0.5;
W = 0.04;
Uc = 1.0e-4;

phiA = 1.0;
phiB = -1.0;

rhoA = 100.0;
rhoB = 1.0;

nuA = 0.1;
nuB = 0.1;

muA = rhoA * nuA;
muB = rhoB * nuB;

fprintf('=============================================\n');
fprintf('Layered Poiseuille analytical check\n');
fprintf('=============================================\n');

fprintf('rhoA/rhoB = %.0f\n', rhoA/rhoB);
fprintf('muA       = %.12e\n', muA);
fprintf('muB       = %.12e\n', muB);
fprintf('Uc        = %.12e\n', Uc);
fprintf('h         = %.12e\n', h);
fprintf('W         = %.12e\n\n', W);


% ------------------------------------------------------------
% 2. Eq.(49): phase field
% ------------------------------------------------------------

phi_fun = @(y) 0.5*(phiA+phiB) + ...
               0.5*(phiA-phiB).*tanh(-2.0.*y./W);


% ------------------------------------------------------------
% 3. Density
% ------------------------------------------------------------

rho_fun = @(y) (rhoA-rhoB)/(phiA-phiB) .* ...
               (phi_fun(y)-phiB) + rhoB;


% ------------------------------------------------------------
% 4. Dynamic viscosity
% ------------------------------------------------------------

mu_fun = @(y) (muA-muB)/(phiA-phiB) .* ...
              (phi_fun(y)-phiB) + muB;


% ------------------------------------------------------------
% 5. Eq.(51): A, B, C, D
% ------------------------------------------------------------

A = integral(@(y) y./mu_fun(y), -h, h, ...
    'AbsTol',1e-13,'RelTol',1e-12);

B = integral(@(y) 1.0./mu_fun(y), -h, h, ...
    'AbsTol',1e-13,'RelTol',1e-12);

C = integral(@(y) y./mu_fun(y), -h, 0, ...
    'AbsTol',1e-13,'RelTol',1e-12);

D = integral(@(y) 1.0./mu_fun(y), -h, 0, ...
    'AbsTol',1e-13,'RelTol',1e-12);

fprintf('A = %.15e\n', A);
fprintf('B = %.15e\n', B);
fprintf('C = %.15e\n', C);
fprintf('D = %.15e\n\n', D);


% ------------------------------------------------------------
% 6. Body force Gx
% ------------------------------------------------------------

Gx = Uc / (A*D/B - C);

fprintf('Gx = %.15e\n\n', Gx);


% ------------------------------------------------------------
% 7. Verify the two coefficients in Eq.(50a)
% ------------------------------------------------------------

coef2_1 = Gx * A / B;
coef2_2 = Uc / (D - B*C/A);

fprintf('Second coefficient check\n');
fprintf('Gx*A/B       = %.15e\n', coef2_1);
fprintf('Paper form   = %.15e\n', coef2_2);
fprintf('difference   = %.15e\n\n', abs(coef2_1-coef2_2));


% ------------------------------------------------------------
% 8. Cell-centered computational grid
% ------------------------------------------------------------

Ny = 100;
dy = 2*h/Ny;

yCell = zeros(Ny,1);

for j = 1:Ny
    yCell(j) = -h + (j-0.5)*dy;
end


% ------------------------------------------------------------
% 9. Analytical u and du/dy
% ------------------------------------------------------------

uExact = zeros(Ny,1);
duExact = zeros(Ny,1);

for j = 1:Ny

    yj = yCell(j);

    I1 = integral(@(yy) yy./mu_fun(yy), -h, yj, ...
        'AbsTol',1e-13,'RelTol',1e-12);

    I2 = integral(@(yy) 1.0./mu_fun(yy), -h, yj, ...
        'AbsTol',1e-13,'RelTol',1e-12);

    % Eq.(50a)
    uExact(j) = -Gx*I1 + Gx*A/B*I2;

    % Eq.(50b)
    duExact(j) = -Gx*yj/mu_fun(yj) + ...
                  Gx*A/B/mu_fun(yj);

end


% ------------------------------------------------------------
% 10. Check u(0) = Uc
% ------------------------------------------------------------

I1_center = integral(@(yy) yy./mu_fun(yy), -h, 0, ...
    'AbsTol',1e-13,'RelTol',1e-12);

I2_center = integral(@(yy) 1.0./mu_fun(yy), -h, 0, ...
    'AbsTol',1e-13,'RelTol',1e-12);

uCenter = -Gx*I1_center + Gx*A/B*I2_center;


% ------------------------------------------------------------
% 11. Check walls
% ------------------------------------------------------------

uBottom = 0.0;

I1_top = integral(@(yy) yy./mu_fun(yy), -h, h, ...
    'AbsTol',1e-13,'RelTol',1e-12);

I2_top = integral(@(yy) 1.0./mu_fun(yy), -h, h, ...
    'AbsTol',1e-13,'RelTol',1e-12);

uTop = -Gx*I1_top + Gx*A/B*I2_top;


fprintf('=============================================\n');
fprintf('Important checks\n');
fprintf('=============================================\n');

fprintf('u_exact(0)  = %.15e\n', uCenter);
fprintf('target Uc   = %.15e\n', Uc);
fprintf('|u(0)-Uc|   = %.15e\n', abs(uCenter-Uc));

fprintf('u(-h)       = %.15e\n', uBottom);
fprintf('u(+h)       = %.15e\n', uTop);

[maxUcell,indexMax] = max(uExact);

fprintf('max cell u  = %.15e\n', maxUcell);
fprintf('y at max u  = %.15e\n', yCell(indexMax));

fprintf('=============================================\n');


% ------------------------------------------------------------
% 12. Fine grid for plotting
% ------------------------------------------------------------

yFine = linspace(-h,h,2001);


% ------------------------------------------------------------
% 13. Plot phi
% ------------------------------------------------------------

figure;

plot(yFine,phi_fun(yFine),'LineWidth',1.6);

xlabel('y');
ylabel('\phi');

title('Initial order parameter');

grid on;
box on;


% ------------------------------------------------------------
% 14. Plot dynamic viscosity
% ------------------------------------------------------------

figure;

plot(yFine,mu_fun(yFine),'LineWidth',1.6);

xlabel('y');
ylabel('\mu(y)');

title('Dynamic viscosity');

grid on;
box on;


% ------------------------------------------------------------
% 15. Plot analytical velocity
% ------------------------------------------------------------

figure;

plot(uExact,yCell,'LineWidth',1.8);

xlabel('u_e');
ylabel('y');

title('Analytical layered Poiseuille velocity');

grid on;
box on;


% ------------------------------------------------------------
% 16. Plot analytical velocity gradient
% ------------------------------------------------------------

figure;

plot(duExact,yCell,'LineWidth',1.8);

xlabel('\partial u_e/\partial y');
ylabel('y');

title('Analytical velocity gradient');

grid on;
box on;


% ------------------------------------------------------------
% 17. Save analytical solution
%
% columns:
% y  phi  rho  mu  uExact  duExact
% ------------------------------------------------------------

phiCell = phi_fun(yCell);
rhoCell = rho_fun(yCell);
muCell = mu_fun(yCell);

reference = [ ...
    yCell, ...
    phiCell, ...
    rhoCell, ...
    muCell, ...
    uExact, ...
    duExact ...
    ];

writematrix(reference, ...
    'poiseuille_exact_100to1.dat', ...
    'Delimiter','space');

fprintf('\nSaved: poiseuille_exact_100to1.dat\n');