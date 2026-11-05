clear;
clc;
close all;

NX = 150;
NY = 600;
d  = 150;

files = {
    'RTI_phi_ACE_Re90_step_05477.dat'
    'RTI_phi_ACE_Re90_step_09585.dat'
    'RTI_phi_ACE_Re90_step_13693.dat'
    'RTI_phi_ACE_Re90_step_16432.dat'
    };

Tlist = [1.0, 1.75, 2.5, 3.0];

x = linspace(-d/2 + 0.5, d/2 - 0.5, NX);
y = linspace(-2*d + 0.5, 2*d - 0.5, NY);

figure('Position',[100 100 1200 700]);

for k = 1:4

    phi = load(files{k});

    subplot(1,4,k);

    imagesc(x,y,phi);

    set(gca,'YDir','normal');

    axis equal;
    axis tight;

    hold on;

    contour(x, y, phi,[0.5 0.5],'k','LineWidth', 1.5);

    caxis([0 1]);

    xlabel('x');

    if k == 1
        ylabel('y');
    end

    title(sprintf('T = %.2f',Tlist(k)));

end

sgtitle('2D RTI, ACE, Re = 90');

colormap(parula);
colorbar('Position',[0.92 0.18 0.015 0.68]);