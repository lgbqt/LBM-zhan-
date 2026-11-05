clear;
clc;
close all;

%==============================
% 你目前实际拥有的文件
%==============================
steps = [0 100000 200000 300000 400000 500000];

colors = lines(length(steps));

figure;
hold on;

for n = 1:length(steps)

    step = steps(n);

    % 文件名
    % step = 0      -> phi_00000.dat
    % step = 100000 -> phi_100000.dat
    filename = sprintf('phi_%05d.dat', step);

    if ~isfile(filename)
        warning('文件 %s 不存在，跳过。', filename);
        continue;
    end

    % 读取 200×200 的 phi 矩阵
    phi = load(filename);

    %==============================
    % 建立物理坐标
    %==============================
    [Ny, Nx] = size(phi);

    dx = 2.0 / Nx;
    dy = 2.0 / Ny;

    x = linspace(-1 + dx/2, 1 - dx/2, Nx);
    y = linspace(-1 + dy/2, 1 - dy/2, Ny);

    %==============================
    % 画 phi = 0 界面
    %==============================
    contour(x, y, phi, [0 0], ...
        'LineWidth', 1.8, ...
        'LineColor', colors(n,:));

end

%==============================
% 图形设置
%==============================
axis equal;

xlim([-0.65 0.65]);
ylim([-0.65 0.65]);

xlabel('x');
ylabel('y');

title('Square droplet evolution (\phi = 0 contour)');

legend( ...
    't = 0', ...
    't = 10', ...
    't = 20', ...
    't = 30', ...
    't = 40', ...
    't = 50', ...
    'Location', 'eastoutside');

grid on;
box on;

hold off;