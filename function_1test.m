k = 10;
x_values = linspace(0, k, 1000);
y_values = function_1(x_values, k);

figure;
plot(x_values, y_values, 'LineWidth', 2);
title('Custom Function');
xlabel('x');
ylabel('f(x)');
grid on;