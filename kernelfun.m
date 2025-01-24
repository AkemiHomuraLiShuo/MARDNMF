% 计算在根号二分之一处的函数值
result_at_sqrt_half = sqrt(2)/2 * exp(-(sqrt(2)/2)^2);

% 计算在1.7338处的函数值
result_at_1_7338 = 1.7338 * exp(-(1.7338)^2);

% 显示结果
disp(['在根号二分之一处的函数值：', num2str(result_at_sqrt_half)]);
disp(['在1.7338处的函数值：', num2str(result_at_1_7338)]);

