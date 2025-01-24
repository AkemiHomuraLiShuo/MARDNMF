% g_function.m
function result = function_2(x, N, k_2)
    condition1 = 0 <= x;
    condition2 = x < N - 1/3 * k_2;
    
    result = zeros(size(x)); % 初始化与 x 相同大小的结果向量
    
    idx = condition1 & condition2;
    result(idx) = 1 - 1/3 * (2.^(3/k_2 * (x(idx) - N + k_2)) - 1);
    
    idx = N - 1/3 * k_2 <= x & x < N;
    result(idx) = 1-1;
    
    % 将结果中大于1的部分设置为1
    result(result > 1) = 1;
end

