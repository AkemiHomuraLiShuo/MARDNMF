function assignment = hungarian_algorithm(matching_matrix)
    % 匈牙利算法实现
    % 输入参数 matching_matrix 是匹配矩阵，其中每个元素表示将一个聚类与一个真实标签匹配的样本数量。

    % 步骤1：减去每行的最小值
    cost_matrix = bsxfun(@minus, matching_matrix, min(matching_matrix, [], 2));

    % 步骤2：减去每列的最小值
    cost_matrix = bsxfun(@minus, cost_matrix, min(cost_matrix, [], 1));

    % 步骤3：找出所有零点
    [row, col] = find(cost_matrix == 0);

    % 步骤4：将零点连接起来
    assignment = zeros(size(cost_matrix));

    % 遍历每一行，为每行找到第一个零点
    for i = 1:size(assignment, 1)
        zero_col = find(row == i, 1);
        if ~isempty(zero_col)
            assignment(i, col(zero_col)) = 1;
        end
    end

    % 遍历每一列，为每列找到第一个零点
    for j = 1:size(assignment, 2)
        zero_row = find(col == j, 1);
        if ~isempty(zero_row)
            assignment(row(zero_row), j) = 1;
        end
    end
end




