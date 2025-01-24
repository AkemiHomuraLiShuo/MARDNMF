%读取数据
load('D:\matlab\donnotquit\dataset\3-sources.mat')
%将每个标签对应数据的位置找出来
indices_1 = find(truth == 1);
indices_2 = find(truth == 2);
indices_3 = find(truth == 3);
indices_4 = find(truth == 4);
indices_5 = find(truth == 5);
indices_6 = find(truth == 6);
%每个簇中只选择其中前11个样本的位置
new_vector_1 = indices_1(1:11);
new_vector_2 = indices_2(1:11);
new_vector_3 = indices_3(1:11);
new_vector_4 = indices_4(1:11);
new_vector_5 = indices_5(1:11);
new_vector_6 = indices_6(1:11);
%将bbc中每个簇对应的元素提取出来
new_bbc1=bbc(new_vector_1, :);
new_bbc2=bbc(new_vector_2, :);
new_bbc3=bbc(new_vector_3, :);
new_bbc4=bbc(new_vector_4, :);
new_bbc5=bbc(new_vector_5, :);
new_bbc6=bbc(new_vector_6, :);
%将这6个数据集中的内容合并为一个数据集
reduced_bbc = vertcat(new_bbc1, new_bbc2, new_bbc3, new_bbc4,new_bbc5,new_bbc6);

%将guardian中每个簇对应的元素提取出来
new_guardian1=guardian(new_vector_1, :);
new_guardian2=guardian(new_vector_2, :);
new_guardian3=guardian(new_vector_3, :);
new_guardian4=guardian(new_vector_4, :);
new_guardian5=guardian(new_vector_5, :);
new_guardian6=guardian(new_vector_6, :);
%将这6个数据集中的内容合并为一个数据集
reduced_guardian = vertcat(new_guardian1, new_guardian2, new_guardian3, new_guardian4,new_guardian5,new_guardian6);

new_reuters1=reuters(new_vector_1, :);
new_reuters2=reuters(new_vector_2, :);
new_reuters3=reuters(new_vector_3, :);
new_reuters4=reuters(new_vector_4, :);
new_reuters5=reuters(new_vector_5, :);
new_reuters6=reuters(new_vector_6, :);
%将这6个数据集中的内容合并为一个数据集
reduced_reuters = vertcat(new_reuters1, new_reuters2, new_reuters3, new_reuters4,new_reuters5,new_reuters6);

%%这个数据集的正确标签是一个列标签
% 定义数字序列和重复次数
numbers = [1, 2, 3, 4, 5, 6];
repetitions = 11;

% 生成列向量
new_truth = reshape(repmat(numbers, repetitions, 1), [], 1);




