% 设定参数
N = 66;
k2 = 33;
k = 10;
% 生成 x 的值
x_values = linspace(0, N-0.001, 1000);
x_values2= linspace(0, k, 1000);
% 计算对应的 y 值
y_values = function3(x_values, N, k2)
y_values2 = 0.5*function_1(x_values2,k)
% 绘制图形
figure;
plot(x_values, y_values, 'LineWidth', 2);
hold on;
plot(x_values2, y_values2, 'LineWidth', 2);
%title('Custom Function g(x)');
set(gca, 'FontSize', 16, 'LineWidth', 2);
xlim([0 66]);
xlabel('x');
ylabel('g(x)');
grid on;
