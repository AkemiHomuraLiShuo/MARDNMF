%1.对数据进行一次kmeans聚类

% 假设 data 是你的矩阵，每一行是一个样本
k = 6;  % 聚类的数目
% 执行 K-Means 聚类
[idx1, centers1] = kmeans(m_scaled_bbc, k);
[idx2, centers2] = kmeans(m_scaled_bbc, k);
[idx3, centers3] = kmeans(m_scaled_bbc, k);

% idx 是一个列向量，表示每个样本所属的聚类
% centers 是一个矩阵，每一行是一个聚类中心的坐标
%1.2计算聚类精度
% 假设 new_truth 是你的真实标签，idx 是聚类结果
%构造函数new_truth
%aidx=idx
best_idx1 = map_clusters(idx1, new_truth, k);
best_idx2 = map_clusters(idx2, new_truth, k);
best_idx3 = map_clusters(idx3, new_truth, k);
%best_idx的作用是将聚类标签尽可能贴近真实情况

%展示上述匹配过程中的matching_matrix
 matching_matrix1 = zeros(k, k);
    for i = 1:k
        for j = 1:k
            matching_matrix1(i, j) = sum((idx1 == i) & (new_truth == j));
        end
    end
 matching_matrix2 = zeros(k, k);
    for i = 1:k
        for j = 1:k
            matching_matrix2(i, j) = sum((idx2 == i) & (new_truth == j));
        end
    end
 matching_matrix3 = zeros(k, k);
    for i = 1:k
        for j = 1:k
            matching_matrix3(i, j) = sum((idx3 == i) & (new_truth == j));
        end
    end    
    
    
    
% 计算 ACC
correct_matches1 = sum(best_idx1 == new_truth);
total_samples = length(new_truth);
ACC1 = correct_matches1 / total_samples;

correct_matches2 = sum(best_idx2 == new_truth);
total_samples = length(new_truth);
ACC2 = correct_matches2 / total_samples;

correct_matches3 = sum(best_idx3 == new_truth);
total_samples = length(new_truth);
ACC3 = correct_matches3 / total_samples;

% 显示 ACC
disp(['Accuracy (ACC1): ', num2str(ACC1)]);
disp(['Accuracy (ACC2): ', num2str(ACC2)]);
disp(['Accuracy (ACC3): ', num2str(ACC3)]);
% 基于best_idx 计算矩阵




%标题计算F矩阵
%首先初始化PreF矩阵的尺度
PreF1 = zeros(length(best_idx1), k);
PreF2 = zeros(length(best_idx2), k);
PreF3 = zeros(length(best_idx3), k);
% 将对应位置设置为 1
for i = 1:length(best_idx1)
    PreF1(i, best_idx1(i)) = 1;
end
for i = 1:length(best_idx2)
    PreF2(i, best_idx2(i)) = 1;
end
for i = 1:length(best_idx3)
    PreF3(i, best_idx3(i)) = 1;
end



%标题，计算PreB矩阵
% 获取每一列非零元素的位置并存储为列向量，这里先处理bbc的任务
num_cols1 = size(PreF1, 2);
nonzero_positions_cell1 = cell(1, num_cols1);%这是一个cell数据，每行是6个元素，每个元素中储存了一些指标
for col = 1:num_cols1
    [row_indices, ~, ~] = find(PreF1(:, col));
    nonzero_positions_cell1{col} = row_indices;
end
% 将列向量按顺序排列成矩阵 PreB
max_nonzero_count1 = max(cellfun(@numel, nonzero_positions_cell1));%这个函数似乎找到了最大容量的Cell中的元素
PreB1 = zeros(max_nonzero_count1, num_cols1);%依然使用之前的列
for col = 1:num_cols1
    PreB1(1:numel(nonzero_positions_cell1{col}), col) = nonzero_positions_cell1{col};%对PreB的每一列进行赋值
end

%重复上述过程处理reuters
num_cols2 = size(PreF2, 2);
nonzero_positions_cell2 = cell(1, num_cols2);%这是一个cell数据，每行是6个元素，每个元素中储存了一些指标
for col = 1:num_cols2
    [row_indices, ~, ~] = find(PreF2(:, col));
    nonzero_positions_cell2{col} = row_indices;
end
% 将列向量按顺序排列成矩阵 PreB2
max_nonzero_count2 = max(cellfun(@numel, nonzero_positions_cell2));%这个函数似乎找到了最大容量的Cell中的元素
PreB2 = zeros(max_nonzero_count2, num_cols2);%依然使用之前的列
for col = 1:num_cols2
    PreB2(1:numel(nonzero_positions_cell2{col}), col) = nonzero_positions_cell2{col};%对PreB的每一列进行赋值
end
%重复上述过程处理guardian
num_cols3 = size(PreF3, 2);
nonzero_positions_cell3 = cell(1, num_cols3);%这是一个cell数据，每行是6个元素，每个元素中储存了一些指标
for col = 1:num_cols3
    [row_indices, ~, ~] = find(PreF3(:, col));
    nonzero_positions_cell3{col} = row_indices;
end
% 将列向量按顺序排列成矩阵 PreB3
max_nonzero_count3 = max(cellfun(@numel, nonzero_positions_cell3));%这个函数似乎找到了最大容量的Cell中的元素
PreB3 = zeros(max_nonzero_count3, num_cols3);%依然使用之前的列
for col = 1:num_cols3
    PreB3(1:numel(nonzero_positions_cell3{col}), col) = nonzero_positions_cell3{col};%对PreB的每一列进行赋值
end



%大标题，计算B矩阵
%接下来我们需要计算出初始化向量的中心preBB一共有6列，每列中都有很多数值，我们需要根据这些数值计算中心向量
%对reduced_bbc进行处理并进行计算得到首个初始化的B
% 获取矩阵的大小
[~, cols_PreB1] = size(PreB1);
% 初始化结果矩阵 Bbb
B1 = zeros(size(reduced_bbc, 2), size(PreB1, 2));
% 计算选择行的平均向量
for i = 1:cols_PreB1
    row_indices = PreB1(:, i);
    
    % 排除0
    valid_indices = row_indices(row_indices > 0);
    
    if ~isempty(valid_indices)
        selected_rows = reduced_bbc(valid_indices, :);
        average_vector = mean(selected_rows, 1);
        B1(:, i) = average_vector';
    else
        % 在 row_indices 中没有正整数的情况下的处理
        % 这里可以根据你的具体需求来决定怎么处理
        % 例如，可以将 Bbb(:, i) 设置为特定的值或者留空
    end
end

% 获取矩阵的大小
[~, cols_PreB2] = size(PreB2);
% 初始化结果矩阵 Bbb
B2 = zeros(size(reduced_reuters, 2), size(PreB2, 2));
% 计算选择行的平均向量
for i = 1:cols_PreB2
    row_indices = PreB2(:, i);
    
    % 排除0
    valid_indices = row_indices(row_indices > 0);
    
    if ~isempty(valid_indices)
        selected_rows = reduced_reuters(valid_indices, :);
        average_vector = mean(selected_rows, 1);
        B2(:, i) = average_vector';
    else
        % 在 row_indices 中没有正整数的情况下的处理
        % 这里可以根据你的具体需求来决定怎么处理
        % 例如，可以将 Bbb(:, i) 设置为特定的值或者留空
    end
end

% 获取矩阵的大小
[~, cols_PreB3] = size(PreB3);
% 初始化结果矩阵 Bbb
B3 = zeros(size(reduced_guardian, 2), size(PreB3, 2));
% 计算选择行的平均向量
for i = 1:cols_PreB3
    row_indices = PreB3(:, i);
    
    % 排除0
    valid_indices = row_indices(row_indices > 0);
    
    if ~isempty(valid_indices)
        selected_rows = reduced_guardian(valid_indices, :);
        average_vector = mean(selected_rows, 1);
        B3(:, i) = average_vector';
    else
        % 在 row_indices 中没有正整数的情况下的处理
        % 这里可以根据你的具体需求来决定怎么处理
        % 例如，可以将 Bbb(:, i) 设置为特定的值或者留空
    end
end




%至此我们已经完成了 1.C矩阵的设计 2.B的初始化，以及3.F的初始化
%2. 根据聚类结果来初始化B,F矩阵，在这个初始化的过程中我们要保证初始化的B的列是基的表达形式       完成
%3问题，我们需要一种新的算法来进行多模态聚类这涉及到kmeans的多模态情况。
%3.根据初始化的B,F矩阵来计算C矩阵，Q，矩阵，以及F*矩阵。其中C矩阵需要在第一轮迭代之后再进行



%标题计算Q矩阵以及F*矩阵，要注意再我们的非负矩阵分解过程中
%Q矩阵是每个B矩阵的列和。
%计算B1对应的Q1
Q1=diag(sum(B1));
%计算B2对应的Q2
Q2=diag(sum(B2));
%计算B3对应的Q3
Q3=diag(sum(B3));

%计算规范化后的F，其实就是F*Q,并按照一范数进行加权求和，除此之外还需要一个beta0
beta0=0.1;%beta0的初始化,这个参数是沟通项和其他项之间的权重比值

beta1 = sum(abs(m_scaled_bbc(:)));
beta2 = sum(abs(m_scaled_reuters(:)));
beta3 = sum(abs(m_scaled_guardian(:)));%这三个数值是各个模态之间的比值
%各个视角下的视角平衡参数设置为1
alpha1=1;
alpha2=1;
alpha3=1;
%每个视角下的损失函数和图正则化项的平衡项参数设置为：
alpha0=1;%目前这个参数还尚不清楚应该设置为多少
%或者说应该将alpha1,alpha2,alpha3,这三个参数分别和图正则化项对应起来
%关于超图正则化项方面还有一些问题
Fx = beta1*PreF1*Q1+beta2*PreF2*Q2+beta3*PreF3*Q3  %Fx是协同矩阵
%注意在buildingC中还有两个关键参数反别是 kdim=0.97 以及 gamma=0.4


%现在找到所有需要调节的参数陈列如下
%1,kdim 2,alpha0图正则化项的系数 3,gamma[0,1]两个核函数的混合平衡系数 4,beta0[0.1,100]沟通项



