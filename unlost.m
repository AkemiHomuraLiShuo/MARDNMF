k = 6;  % 聚类的数目
% 构建匹配矩阵，计算每个聚类与真实标签的匹配成本
matching_matrix = zeros(k, k);
for i = 1:k
    for j = 1:k
        matching_matrix(i, j) = sum((idx == i) & (new_truth == j));
    end
end
% 使用匈牙利算法找到最优匹配(这里使用的是一个)
G= munkers(-matching_matrix)% 示例最佳匹配矩阵
% 获取分配矩阵的行索引或列索引
[~, column_indices] = max(G, [], 2);
% 使用索引映射原始向量的元素
best_idx = zeros(size(idx));
for i = 1:length(column_indices)
    best_idx(idx == i) = column_indices(i);
end