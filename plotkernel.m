% 定义符号变量 x
syms x

% 定义函数 f(x) = xe^(-x^2)
f = (2/9)*x*2*exp(-x^2/9);
f2 = 0.125*x*4*exp(-(x^2/16));
f3 = f + f2;

% 使用 ezplot 绘制图像
figure;
hold on;  % 保持图形，使得可以在同一图中添加多个曲线

% 绘制 f(x)
h1 = ezplot(f, [0, 16]);
set(h1, 'LineWidth', 2);  % 设置线条宽度

% 绘制 f2(x)
h2 = ezplot(f2, [0, 16]);
set(h2, 'LineWidth', 2);  % 设置线条宽度

% 绘制 f3(x)
h3 = ezplot(f3, [0, 16]);
set(h3, 'LineWidth', 2);  % 设置线条宽度

hold off;  % 取消图形保持
set(gcf, 'NumberTitle', 'off');

% 移除图像标题
title('');

% 添加标签
xlabel('x', 'FontSize', 16);  % 设置 x 轴标签的字体大小
ylabel('y', 'FontSize', 16);  % 设置 y 轴标签的字体大小

% 设置图例
legend([h1, h2, h3], {'f(x)', 'f2(x)', 'f3(x)'}, 'Location', 'Best', 'FontSize', 15, 'FontWeight', 'normal');  % 设置图例字体大小和粗细


% 设置坐标轴刻度标签和坐标轴的字体大小和宽度
set(gca, 'FontSize', 16, 'LineWidth', 2);  %


