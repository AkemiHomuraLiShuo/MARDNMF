function y = function_1(x, k)
    y = zeros(size(x));

    % 第一个分段
    idx1 = (0 <= x) & (x < (2/5)*k);
    y(idx1) =2* 1;

    % 第二个分段
    idx2 = ((2/5)*k <= x) & (x < k);
    y(idx2) = 2* 1/7 * (2.^(5/k * (-x(idx2) + k)) - 1);
end
