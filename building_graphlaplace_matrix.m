load('reduced3sources.mat')
% 对数据进行缩放，按照二范数进行统一

% 计算每个模态每行的二范数
row_norms_bbc = vecnorm(reduced_bbc, 2, 2);
row_norms_reuters = vecnorm(reduced_reuters, 2, 2);
row_norms_guardian = vecnorm(reduced_guardian, 2, 2);
% 计算这里每个模态中样本的平均值
average_value_bbc = mean(row_norms_bbc );
average_value_reuters = mean(row_norms_reuters );
average_value_guardian = mean(row_norms_guardian );
% 计算缩放比例
scale_factors_bbc = average_value_bbc ./ row_norms_bbc;
scale_factors_reuters = average_value_reuters ./ row_norms_reuters;
scale_factors_guardian = average_value_guardian ./ row_norms_guardian;
% 缩放矩阵的每一行使得每行的二范数等于 average_value
m_scaled_bbc = reduced_bbc .* scale_factors_bbc;
m_scaled_reuters = reduced_reuters .* scale_factors_reuters;
m_scaled_guardian = reduced_guardian .* scale_factors_guardian;

%  1  计算三个模态的距离矩阵D，既保留按照二范数缩放过后的，也保留没有按照二范数缩放的。
m_distance_scaledbbc= pdist2(m_scaled_bbc, m_scaled_bbc);
m_distance_scaledreuters= pdist2(m_scaled_reuters, m_scaled_reuters);
m_distance_scaledguardian= pdist2(m_scaled_guardian, m_scaled_guardian);
%没有按照二范数缩放的矩阵如下：

%注意此处应该补充另外三个距离矩阵，针对没有缩放过的数据集
%计算reduced_bbc这个矩阵中样本的距离矩阵D_bbc
%计算reduced_reuters这个矩阵中样本的距离矩阵D_reuters
%计算reduced_guardian这个矩阵中样本的距离矩阵D_guardian




%2.   计算由距离矩阵得到的次序矩阵R，
%    2.1 计算bbc的次序矩阵
[~, sorted_indices] = sort(m_distance_scaledbbc);
for col = 1:size(m_distance_scaledbbc, 2)
    [~, order] = sort(sorted_indices(:, col));
    MRscaledbbc(:, col) = order;
end
%    2.2 计算reuters的次序矩阵
[~, sorted_indices] = sort(m_distance_scaledreuters);
for col = 1:size(m_distance_scaledreuters, 2)
    [~, order] = sort(sorted_indices(:, col));
    MRscaledreuters(:, col) = order;
end
%    2.3 计算guardian的次序矩阵
[~, sorted_indices] = sort(m_distance_scaledguardian);
for col = 1:size(m_distance_scaledguardian, 2)
    [~, order] = sort(sorted_indices(:, col));
    MRscaledguardian(:, col) = order;
end
%%%%%%%%%      非常重要的一点  设定k值的大小，在运行前需要设置的大小
k=9; k2=33; N=66 ;
%构造函数。函数体名称为function_1和function_2

%  3.   将次序矩阵中数值小于k的样本记为1，同时令对角线元素为0得到指标矩阵H
% 将小于 k 的元素设置为 1，其余设置为 0
MHbbc= zeros(size(MRscaledbbc)); 
MHbbc(MRscaledbbc < k) = 1;
MHreuters= zeros(size(MRscaledreuters)); 
MHreuters(MRscaledreuters < k) = 1;
MHguardian= zeros(size(MRscaledguardian));  
MHguardian(MRscaledguardian < k) = 1;
MHbbc(logical(eye(size(MHbbc)))) = 0;
MHreuters(logical(eye(size(MHreuters)))) = 0;
MHguardian(logical(eye(size(MHguardian)))) = 0;

%将大于k2的元素设置为1，其余设置为0

MH2bbc= zeros(size(MRscaledbbc)); 
MH2bbc(MRscaledbbc > N-k2) = 1;
MH2reuters= zeros(size(MRscaledreuters)); 
MH2reuters(MRscaledreuters > N-k2) = 1;
MH2guardian= zeros(size(MRscaledguardian));  
MH2guardian(MRscaledguardian > N-k2) = 1;
MH2bbc(logical(eye(size(MH2bbc)))) = 0;
MH2reuters(logical(eye(size(MH2reuters)))) = 0;
MH2guardian(logical(eye(size(MH2guardian)))) = 0;




%绘制热图
% 创建热图
heatmap(MHbbc, 'Colormap', hot, 'ColorLimits', [min(MHbbc(:)), max(MHbbc(:))]);
% 添加标题和标签
title('Heatmap of H');
xlabel('Columns');
ylabel('Rows');
% 创建热图
heatmap(MHreuters, 'Colormap', hot, 'ColorLimits', [min(MHreuters(:)), max(MHreuters(:))]);
% 添加标题和标签
title('Heatmap of H');
xlabel('Columns');
ylabel('Rows');
% 创建热图
heatmap(MHguardian, 'Colormap', hot, 'ColorLimits', [min(MHguardian(:)), max(MHguardian(:))]);
% 添加标题和标签
title('Heatmap of H');
xlabel('Columns');
ylabel('Rows');

%  5.   将正向鼓励函数作用次序矩阵R得到一个A矩阵
N=66;K2=33;
MAbbc = function_1(MRscaledbbc, k);
MAreuters = function_1(MRscaledreuters, k);
MAguardian= function_1(MRscaledguardian, k);


% 5.5 将负向惩罚函数作用在次序矩阵R上 得到A2矩阵

MA2bbc = function_2(MRscaledbbc, N, k2);
MA2reuters = function_2(MRscaledreuters, N, k2);
MA2guardian= function_2(MRscaledguardian, N, k2);
heatmap(MA2bbc, 'Colormap', hot, 'ColorLimits', [min(MA2bbc(:)), max(MA2bbc(:))]);
Test1=MA2bbc+MA2bbc';
Test1(Test1 > 2) = 2;
h2=heatmap(Test1, 'Colormap', hot, 'ColorLimits', [min(Test1(:)), max(Test1(:))]);
%MA2的图效果可能更好一些
h = gca;
h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
h.YDisplayLabels = repmat({''}, size(h.YDisplayData));
h2.FontSize=14;
% 隐藏网格线
%h.GridVisible = 'off';

% 隐藏坐标轴
h.XLabel = '';
h.YLabel = '';




%  6.   将A矩阵于H矩阵元素相乘得到一个对角矩阵记为B矩阵
MBbbc=MAbbc .* MHbbc;
MBreuters=MAreuters .* MHreuters;
MBguardian=MAguardian .* MHguardian;
%对角线归零
MBbbc(logical(eye(size(MBbbc)))) = 0;
MBreuters(logical(eye(size(MBreuters)))) = 0;
MBguardian(logical(eye(size(MBguardian)))) = 0;
% 创建热图,整体预览
%heatmap(MBbbc, 'Colormap', hot, 'ColorLimits', [min(MBbbc(:)), max(MBbbc(:))]);
Test2=0.5*(MBbbc+MBbbc');
h1=heatmap(Test2, 'Colormap', hot, 'ColorLimits', [min(Test2(:)), max(Test2(:))]);
% 添加标题和标签
% 去掉X轴、Y轴的标签
h = gca;
h.XDisplayLabels = repmat({''}, size(h.XDisplayData));
h.YDisplayLabels = repmat({''}, size(h.YDisplayData));
h1.FontSize=14;
% 隐藏网格线
%h.GridVisible = 'off';

% 隐藏坐标轴
h.XLabel = '';
h.YLabel = '';




% 权重矩阵Wbbcs Wreuterss Wguardians
%MB也是可以用的
%f1函数和f2函数的性质？
%f1不变，但是f2数值上增加了一个1

%到6.1为止我们重新构建图laplace矩阵，6.1是对鼓励图进行处理，因此是可行的

newlaplacebb=diag(sum(MBbbc))-MBbbc;
newlaplacere=diag(sum(MBreuters))-MBreuters;
newlaplacegu=diag(sum(MBguardian))-MBguardian;

%类似地，我们绘制新的laplace矩阵在第二尺度。基于MA2

MA2bbc(logical(eye(size(MA2bbc)))) = 0;
MA2reuters(logical(eye(size(MA2reuters)))) = 0;
MA2guardian(logical(eye(size(MA2guardian)))) = 0;

newlaplace2bb=diag(sum(MA2bbc))-MA2bbc;
newlaplace2re=diag(sum(MA2reuters))-MA2reuters;
newlaplace2gu=diag(sum(MA2guardian))-MA2guardian;

%%停止






%6.5
MB2bbc=MA2bbc .* MH2bbc;
MB2reuters=MA2reuters .* MH2reuters;
MB2guardian=MA2guardian .* MH2guardian;
%对角线归零
MB2bbc(logical(eye(size(MB2bbc)))) = 0;
MB2reuters(logical(eye(size(MB2reuters)))) = 0;
MB2guardian(logical(eye(size(MB2guardian)))) = 0;
% 创建热图,整体预览
heatmap(MB2bbc, 'Colormap', hot, 'ColorLimits', [min(MB2bbc(:)), max(MB2bbc(:))]);
% 添加标题和标签
title('Heatmap of H');
xlabel('Columns');
ylabel('Rows');
% 权重矩阵Wbbcs Wreuterss Wguardians




%  7.   将B矩阵中的元素计算列和，得到一个对角矩阵即为W矩阵，其中W矩阵对角线上的元素对应着每条边的边权
%计算矩阵 bbc 的列和
col_sums_bbc = sum(MBbbc);
% 创建对角矩阵 W，对角线上的元素为列和
Wbbc = diag(col_sums_bbc);
%计算矩阵 reuters 的列和
col_sums_reuters = sum(MBreuters);
% 创建对角矩阵 W，对角线上的元素为列和
Wreuters = diag(col_sums_reuters);
%计算矩阵 guardian 的列和                                   %出错！！！！
col_sums_guardian = sum(MBguardian);
% 创建对角矩阵 W，对角线上的元素为列和
Wguardian = diag(col_sums_guardian);


% 7.5  重复第7步里的内容我们可以得出惩罚图的设计
%计算矩阵 bbc 的列和
col_sums_2bbc = sum(MB2bbc);
% 创建对角矩阵 W，对角线上的元素为列和
W2bbc = diag(col_sums_2bbc);
%计算矩阵 reuters 的列和
col_sums_2reuters = sum(MB2reuters);
% 创建对角矩阵 W，对角线上的元素为列和
W2reuters = diag(col_sums_2reuters);
%计算矩阵 guardian 的列和                                   %出错！！！！
col_sums_2guardian = sum(MB2guardian);
% 创建对角矩阵 W，对角线上的元素为列和
W2guardian = diag(col_sums_2guardian);








%  8.   将H矩阵与W矩阵相乘，再计算行求和即得到最终的D_H矩阵
DHbbc= diag(sum( MHbbc * Wbbc, 2));
DHreuters= diag(sum(MHreuters * Wreuters, 2));
DHguardian= diag(sum( MHguardian * Wguardian , 2));

%  8.5  负图构建
DH2bbc= diag(sum( MH2bbc * W2bbc, 2));
DH2reuters= diag(sum(MH2reuters * W2reuters, 2));
DH2guardian= diag(sum( MH2guardian * W2guardian , 2));

%  9.   将H矩阵计算列和得到D_e矩阵，其中元素的边的度                  
DEbbc = diag(sum(MHbbc));
DEreuters = diag(sum(MHreuters));
DEguardian = diag(sum(MHguardian));

%  9.5  9的负图版本
DE2bbc = diag(sum(MH2bbc));
DE2reuters = diag(sum(MH2reuters));
DE2guardian = diag(sum(MH2guardian));


%  10.  最后得到图laplace矩阵通过D_H-HWD_e^{-1}H^T=D_H-S_H
Laplacianbbc=DHbbc - MHbbc * Wbbc * inv(DEbbc)* MHbbc';
sum(Laplacianbbc);
Laplacianreuters=DHreuters - MHreuters * Wreuters * inv(DEreuters)* MHreuters';
sum(Laplacianreuters);
Laplacianguardian=DHguardian - MHguardian * Wguardian * inv(DEguardian)* MHguardian';
sum(Laplacianguardian);
%218行

%  10.5  最后得到图laplace矩阵通过D_H-HWD_e^{-1}H^T=D_H-S_H
Laplacian2bbc=DH2bbc - MH2bbc * W2bbc * inv(DE2bbc)* MH2bbc';
a=sum(Laplacian2bbc);
Laplacian2reuters=DH2reuters - MH2reuters * W2reuters * inv(DE2reuters)* MH2reuters';
b=sum(Laplacian2reuters);
Laplacian2guardian=DH2guardian - MH2guardian * W2guardian * inv(DE2guardian)* MH2guardian';
c=sum(Laplacian2guardian);


%modified_matrix = Laplacianbbc - diag(diag(Laplacianbbc));
%heatmap(modified_matrix, 'Colormap', hot, 'ColorLimits', [min(modified_matrix(:)), max(modified_matrix(:))]);



