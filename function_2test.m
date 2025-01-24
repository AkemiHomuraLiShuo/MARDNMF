% 设定参数
N = 66;
k2 = 33;

% 生成 x 的值
x_values = linspace(1, N-0.001, 1000);

% 计算对应的 y 值
y_values = function_2(x_values, N, k2)

% 绘制图形
figure;
plot(x_values, y_values, 'LineWidth', 2);
hold on;
title('Custom Function g(x)');
xlabel('x');
ylabel('g(x)');
grid on;

% 调整 y 轴范围
%ylim([-2, 2]); % 适当选择范围
%x=78
%y=function_2(x, N, k2);