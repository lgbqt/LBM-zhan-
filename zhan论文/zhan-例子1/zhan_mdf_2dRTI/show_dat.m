clear;
clc;
close all;

NX = 150;
NY = 600;

d = 150;

filename = 'RTI_phi_ACE_Re90_step_05477.dat';

phi = load(filename);

x = linspace(-d/2 + 0.5, d/2 - 0.5, NX);
y = linspace(-2*d + 0.5, 2*d - 0.5, NY);

figure;

imagesc(x, y, phi);

set(gca, 'YDir', 'normal');

axis equal;
axis tight;

colorbar;

hold on;

contour(x,y,phi,[0.5 0.5],'k','LineWidth',2);

xlabel('x');
ylabel('y');

title('2D RTI, ACE, Re = 90, T \approx 0.365');