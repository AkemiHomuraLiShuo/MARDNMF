% g_function.m
function result = function3(x, N, k_2)
    % 初始化与 x 相同大小的结果向量
    result = zeros(size(x)); 
    
    % 第一段： 0 <= x < N - k_2，结果取值为 1
    condition1 = 0 <= x & x < N - k_2;
    result(condition1) = 1;
    
    % 第二段： N - k_2 <= x < N - 1/3 * k_2，执行指数计算
    condition2 = N - k_2 <= x & x < N - 1/3 * k_2;
    result(condition2) = 1 - 1/3 * (2.^(3/k_2 * (x(condition2) - N + k_2)) - 1);
    
    % 第三段： N - 1/3 * k_2 <= x < N，结果取值为 0
    condition3 = N - 1/3 * k_2 <= x & x < N;
    result(condition3) = 0; % 等价于 1 - 1;
    
    % 其他区间可以根据需要继续扩展
end
