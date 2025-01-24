
%      prework
load('proceed data.mat')




%      1 对三个模态或者视角数据集进行切割，切割成6个部分，每个模态会形成6个子类

%1.1计算bbc的拆分
% 拆分矩阵
submatricesbbc1 = struct();
num_rows_per_submatrix = 11;
num_submatrices = size(m_scaled_bbc, 1) / num_rows_per_submatrix;

% 创建子矩阵并存储在结构体数组中
for i = 1:num_submatrices
    start_row = (i - 1) * num_rows_per_submatrix + 1;
    end_row = i * num_rows_per_submatrix;
    submatrix_name = ['Mbbc', num2str(i)];
    
    % 存储子矩阵在结构体数组中
    submatricesbbc1.(submatrix_name) = m_scaled_bbc(start_row:end_row, :);
end

%计算reuters的拆分
% 拆分矩阵
submatricesreuters2 = struct();
num_rows_per_submatrix = 11;
num_submatrices = size(m_scaled_reuters, 1) / num_rows_per_submatrix;

% 创建子矩阵并存储在结构体数组中
for i = 1:num_submatrices
    start_row = (i - 1) * num_rows_per_submatrix + 1;
    end_row = i * num_rows_per_submatrix;
    submatrix_name = ['Mreuters', num2str(i)];
    
    % 存储子矩阵在结构体数组中
    submatricesreuters2.(submatrix_name) = m_scaled_reuters(start_row:end_row, :);
end

%计算guardian的拆分
% 拆分矩阵
submatricesguardian3 = struct();
num_rows_per_submatrix = 11;
num_submatrices = size(m_scaled_guardian, 1) / num_rows_per_submatrix;

% 创建子矩阵并存储在结构体数组中
for i = 1:num_submatrices
    start_row = (i - 1) * num_rows_per_submatrix + 1;
    end_row = i * num_rows_per_submatrix;
    submatrix_name = ['Mguardian', num2str(i)];
    
    % 存储子矩阵在结构体数组中
    submatricesguardian3.(submatrix_name) =m_scaled_guardian(start_row:end_row, :);
end


%      2 按照一个模态来说，上述操作将数据分为了6个子类，我们需要在本阶段中计算子类中样本的距离矩阵

%计算bbc中的分组距离矩阵
% 创建一个空的结构体数组用于存储距离矩阵
Dsubmatbbc = struct();
% 计算每个子结构体中的样本间的距离矩阵
for i = 1:6
    matrix_name = ['Mbbc', num2str(i)];
    matrix_name2 = ['Dbbc', num2str(i)];
    % 使用 pdist 函数计算距离
    distances = pdist(submatricesbbc1.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsubmatbbc(1).(matrix_name2) = squareform(distances); % 将距离转换为距离矩阵形式
end

%计算reuters中的分组距离矩阵
% 创建一个空的结构体数组用于存储距离矩阵
Dsubmatreuters = struct();
% 计算每个子结构体中的样本间的距离矩阵
for i = 1:6
    matrix_name = ['Mreuters', num2str(i)];
    matrix_name2 = ['Dreuters', num2str(i)];
    % 使用 pdist 函数计算距离
    distances = pdist(submatricesreuters2.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsubmatreuters(1).(matrix_name2) = squareform(distances); % 将距离转换为距离矩阵形式
end

%计算guardian中的分组距离矩阵
% 创建一个空的结构体数组用于存储距离矩阵
Dsubmatguardian = struct();
% 计算每个子结构体中的样本间的距离矩阵
for i = 1:6
    matrix_name = ['Mguardian', num2str(i)];
     matrix_name2 = ['Dguardian', num2str(i)];
    % 使用 pdist 函数计算距离
    distances = pdist(submatricesguardian3.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsubmatguardian(1).(matrix_name2) = squareform(distances); % 将距离转换为距离矩阵形式
end



%      3 得到距离矩阵之后进行行求和并绘制热图，只是一条线，但是能找到极端值。
%3.1计算bbc的情况
Dsumbbc = struct();
for i = 1:6
    matrix_name = ['Dbbc', num2str(i)];
     matrix_name2 = ['Dsumbbc', num2str(i)];
    % 使用 pdist 函数计算距离
    sumdis = sum(Dsubmatbbc.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsumbbc.(matrix_name2) = sumdis; % 将距离转换为距离矩阵形式
end
%3.2计算reuters的情况
Dsumbbc = struct();
for i = 1:6
    matrix_name = ['Dbbc', num2str(i)];
     matrix_name2 = ['Dsumbbc', num2str(i)];
    % 使用 pdist 函数计算距离
    sumdis = sum(Dsubmatbbc.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsumbbc.(matrix_name2) = sumdis; % 将距离转换为距离矩阵形式
end
%3.2计算reuters的情况
Dsumreuters = struct();
for i = 1:6
    matrix_name = ['Dreuters', num2str(i)];
     matrix_name2 = ['Dsumreuters', num2str(i)];
    % 使用 pdist 函数计算距离
    sumdis = sum(Dsubmatreuters.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsumreuters.(matrix_name2) = sumdis; % 将距离转换为距离矩阵形式
end
%3.3计算guardian的情况
Dsumguardian = struct();
for i = 1:6
    matrix_name = ['Dguardian', num2str(i)];
     matrix_name2 = ['Dsumguardian', num2str(i)];
    % 使用 pdist 函数计算距离
    sumdis = sum(Dsubmatguardian.(matrix_name));
    
    % 将距离矩阵存储在结构体数组中的不同字段
    Dsumguardian.(matrix_name2) = sumdis; % 将距离转换为距离矩阵形式
end

%3.4 选择bbc中的异常样本记为子向量为b
% 创建一个空的结构体数组 B
indicatorbbc = struct('inb1', [], 'inb2', [], 'inb3', [], 'inb4', [], 'inb5', [], 'inb6', []);
% 对每个向量字段进行处理
for i = 1:6
    % 获取当前向量字段的数据
    current_data = Dsumbbc.(['Dsumbbc' num2str(i)]);
    
    % 计算当前向量的平均值的1.1倍
    average_value = 1.015 * (mean(current_data));
    
    % 创建新的向量，根据条件将元素设置为1或0
    new_vector = (current_data > average_value);
    
    % 将新向量赋值给结构体数组 B 中对应的字段
    indicatorbbc.(['inb' num2str(i)]) = new_vector;
end

%3.5 选择reuters中的异常样本记为子向量为b
% 创建一个空的结构体数组 B
indicatorreuters = struct('inr1', [], 'inr2', [], 'inr3', [], 'inr4', [], 'inr5', [], 'inr6', []);
% 对每个向量字段进行处理
for i = 1:6
    % 获取当前向量字段的数据
    current_data = Dsumreuters.(['Dsumreuters' num2str(i)]);
    
    % 计算当前向量的平均值的1.1倍
    average_value = 1.015 * (mean(current_data));
    
    % 创建新的向量，根据条件将元素设置为1或0
    new_vector = (current_data > average_value);
    
    % 将新向量赋值给结构体数组 B 中对应的字段
    indicatorreuters.(['inr' num2str(i)]) = new_vector;
end

%3.5 选择guardian中的异常样本记为子向量为b
% 创建一个空的结构体数组 B
indicatorguardian = struct('ing1', [], 'ing2', [], 'ing3', [], 'ing4', [], 'ing5', [], 'ing6', []);
% 对每个向量字段进行处理
for i = 1:6
    % 获取当前向量字段的数据
    current_data = Dsumguardian.(['Dsumguardian' num2str(i)]);
    
    % 计算当前向量的平均值的1.1倍
    average_value = 1.02 * (mean(current_data));
    
    % 创建新的向量，根据条件将元素设置为1或0
    new_vector = (current_data > average_value);
    
    % 将新向量赋值给结构体数组 B 中对应的字段
    indicatorguardian.(['ing' num2str(i)]) = new_vector;
end


%      4 在每一个子类bag1-6，每个视角6个子类中都会拆分出一个异常集合B集合，与一个正常集合A。可以记为MbbcA1-6,MbbcB1-6,A代表正常，B代表异常，编号对应子类数目。

% 4.1计算bbc中的异常样本并进行分离
% 创建一个空的结构体数组用于存储分割后的矩阵
submatbb_split = struct();

% 循环遍历结构体数组中的矩阵
for i = 1:6
    matrix_name = ['Mbbc', num2str(i)];
    indicator_name = ['inb', num2str(i)];
    % 提取当前矩阵和相应的向量
    current_matrix = submatricesbbc1.(matrix_name);
    current_indicator = indicatorbbc.(indicator_name);
    % 使用逻辑索引分割矩阵
    matrix_split1 = current_matrix(current_indicator == 0, :);
    matrix_split2 = current_matrix(current_indicator == 1, :);
    % 存储分割后的矩阵到新的结构体数组中
    submatbb_split.([matrix_name, '1']) = matrix_split1;

    submatbb_split.([matrix_name, '2']) = matrix_split2;
end


% 4.2计算reuters中的异常样本并分离它们。submatrices用 indecator 分离为 sbumat split
submatre_split = struct();
% 循环遍历结构体数组中的矩阵
for i = 1:6
    matrix_name = ['Mreuters', num2str(i)];
    indicator_name = ['inr', num2str(i)];
% 提取当前矩阵和相应的向量
    current_matrix = submatricesreuters2.(matrix_name);
    current_indicator = indicatorreuters.(indicator_name);
% 使用逻辑索引分割矩阵
    matrix_split1 = current_matrix(current_indicator == 0, :);
    matrix_split2 = current_matrix(current_indicator == 1, :);
% 存储分割后的矩阵到新的结构体数组中
    submatre_split.([matrix_name, '1']) = matrix_split1;

    submatre_split.([matrix_name, '2']) = matrix_split2;
end


% 4.3计算guardian中的异常样本并分离它们 
submatgu_split = struct();
% 循环遍历结构体数组中的矩阵
for i = 1:6
    matrix_name = ['Mguardian', num2str(i)];
    indicator_name = ['ing', num2str(i)];
% 提取当前矩阵和相应的向量
    current_matrix = submatricesguardian3.(matrix_name);
    current_indicator = indicatorguardian.(indicator_name);
% 使用逻辑索引分割矩阵
    matrix_split1 = current_matrix(current_indicator == 0, :);
    matrix_split2 = current_matrix(current_indicator == 1, :);
% 存储分割后的矩阵到新的结构体数组中
    submatgu_split.([matrix_name, '1']) = matrix_split1;

    submatgu_split.([matrix_name, '2']) = matrix_split2;
end


%      5 计算MbbcA1-6的6个中心记为cMbbcA1-6。以及其他模态同样处理，共得到18个中心。

% 5.1计算bbc中的正常核心向量，用于估计数值。
% 假设 submatbb_split 是一个包含12个矩阵的结构体数组
% 例如，submatbb_split.Mbbc11 对应 Mbbc11，submatbb_split.Mbbc12 对应 Mbbc12，以此类推
% 指定要计算中心的矩阵编号
target_matrices = [11, 21, 31, 41, 51, 61];
% 创建一个空的结构体数组用于存储每个矩阵中样本的中心
submatbb_cen = struct();
% 循环遍历指定的矩阵编号
for i = 1:numel(target_matrices)
    matrix_number = target_matrices(i);
    matrix_name = ['Mbbc', num2str(matrix_number)];

    % 提取当前矩阵
    current_matrix = submatbb_split.(matrix_name);

    % 计算样本中心
    center_vector = mean(current_matrix, 1);

    % 存储样本中心到新的结构体数组
    submatbb_cen.(['cbbc', num2str(matrix_number)]) = center_vector;
end

%5.2计算reuters中的正常核心向量，用于估计数值。
submatre_cen = struct();
% 循环遍历指定的矩阵编号
for i = 1:numel(target_matrices)
    matrix_number = target_matrices(i);
    matrix_name = ['Mreuters', num2str(matrix_number)];

    % 提取当前矩阵
    current_matrix = submatre_split.(matrix_name);

    % 计算样本中心
    center_vector = mean(current_matrix, 1);

    % 存储样本中心到新的结构体数组
    submatre_cen.(['creuters', num2str(matrix_number)]) = center_vector;
end


%5.3计算reuters中的正常核心向量，用于估计数值。
submatgu_cen = struct();
% 循环遍历指定的矩阵编号
for i = 1:numel(target_matrices)
    matrix_number = target_matrices(i);
    matrix_name = ['Mguardian', num2str(matrix_number)];

    % 提取当前矩阵
    current_matrix = submatgu_split.(matrix_name);

    % 计算样本中心
    center_vector = mean(current_matrix, 1);

    % 存储样本中心到新的结构体数组
    submatgu_cen.(['cguardian', num2str(matrix_number)]) = center_vector;
end


%      6 计算MbbcA1-6到cMbbcA1-6的平均距离，记为bbcpretheta1-6。计算MbbcB1-6到cMbbcA1-6的平均距离，记为bbcpreverta1-6

%6.1计算bbc正常样本的每个样本到质心的距离
% 初始化 discombb 结构体数组
discombb = struct();
% 循环计算每个矩阵到相应向量的距离
for i = 1:6
    matrix_name = ['Mbbc', num2str(i * 10 + 1)];
    vector_name = ['cbbc', num2str(i * 10 + 1)];
    matrix_name2  = ['Mbbc', num2str(i * 10 + 2)];
    matrix_name3  = ['Mbbc', num2str(i)];
    % 获取当前矩阵和向量
    current_matrix = submatbb_split.(matrix_name);
    current_matrix2 = submatbb_split.(matrix_name2);
    current_matrix3 = submatricesbbc1.(matrix_name3);
    current_vector = submatbb_cen.(vector_name);

    % 计算距离并存储结果
    distances1 = sqrt(sum((current_matrix - current_vector).^2, 2));distances1 = full(distances1);
    distances2 = sqrt(sum((current_matrix2 - current_vector).^2, 2));distances2 = full(distances2);
    distances3 = sqrt(sum((current_matrix3 - current_vector).^2, 2));distances3 = full(distances3);
    
    discombb.(['disbba', num2str(i)]) = distances1;
    discombb.(['disbbb', num2str(i)]) = distances2;
    discombb.(['disbbtotal', num2str(i)]) = distances3;
end


%6.2计算reuters正常样本的每个样本到质心的距离
% 初始化 discombb 结构体数组
discomre = struct();
% 循环计算每个矩阵到相应向量的距离
for i = 1:6
    matrix_name = ['Mreuters', num2str(i * 10 + 1)];
    vector_name = ['creuters', num2str(i * 10 + 1)];
    matrix_name2  = ['Mreuters', num2str(i * 10 + 2)];
    matrix_name3  = ['Mreuters', num2str(i)];
    % 获取当前矩阵和向量
    current_matrix = submatre_split.(matrix_name);
    current_matrix2 = submatre_split.(matrix_name2);
    current_matrix3 = submatricesreuters2.(matrix_name3);
    current_vector = submatre_cen.(vector_name);

    % 计算距离并存储结果
    distances1 = sqrt(sum((current_matrix - current_vector).^2, 2));distances1 = full(distances1);
    distances2 = sqrt(sum((current_matrix2 - current_vector).^2, 2));distances2 = full(distances2);
    distances3 = sqrt(sum((current_matrix3 - current_vector).^2, 2));distances3 = full(distances3);
    
    discomre.(['disrea', num2str(i)]) = distances1;
    discomre.(['disreb', num2str(i)]) = distances2;
    discomre.(['disretotal', num2str(i)]) = distances3;
end



%6.3计算guardian正常样本的每个样本到质心的距离
% 初始化 discombb 结构体数组
discomgu = struct();
% 循环计算每个矩阵到相应向量的距离
for i = 1:6
    matrix_name = ['Mguardian', num2str(i * 10 + 1)];
    vector_name = ['cguardian', num2str(i * 10 + 1)];
    matrix_name2  = ['Mguardian', num2str(i * 10 + 2)];
    matrix_name3  = ['Mguardian', num2str(i)];
    % 获取当前矩阵和向量
    current_matrix = submatgu_split.(matrix_name);
    current_matrix2 = submatgu_split.(matrix_name2);
    current_matrix3 = submatricesguardian3.(matrix_name3);
    current_vector = submatgu_cen.(vector_name);

    % 计算距离并存储结果
    distances1 = sqrt(sum((current_matrix - current_vector).^2, 2));distances1 = full(distances1);
    distances2 = sqrt(sum((current_matrix2 - current_vector).^2, 2));distances2 = full(distances2);
    distances3 = sqrt(sum((current_matrix3 - current_vector).^2, 2));distances3 = full(distances3);
    
    discomgu.(['disgua', num2str(i)]) = distances1;
    discomgu.(['disgub', num2str(i)]) = distances2;
    discomgu.(['disgutotal', num2str(i)]) = distances3;
end




%接下来准备数据维度转换，其中应该注意到这个空间的维度是6维的，因为进行了一次范数归一化，数据分布在6维空间空间中分布，那么问题集中在这两个差异较大的维度应该
%如何进行取舍？N-1维/6维
%计算矩阵的秩
% 假设 reduced_bbc 是你的稀疏矩阵
% 使用 svds 计算矩阵的秩
%这里应该是10维起步，数值上应该是0.94以上

%6.4计算这些量的平均值。首先，计算bbc的平均值
% 初始化新向量 prethetabb
prethetabb = [];
prevertabb = [];

% 循环计算每个向量的平均值
for i = 1:6
    % 获取当前向量
    current_vector = discombb.(['disbba', num2str(i)]);  % 假设向量保存在结构体的 vector 字段中
    current_vector2 = discombb.(['disbbb', num2str(i)]);  
    % 计算当前向量的平均值
    mean_vector = mean(current_vector);
    mean_vector2 = mean(current_vector2);
    % 将平均值追加到向量 prethetabb 中
    prethetabb = [prethetabb, mean_vector];
    prevertabb = [prevertabb, mean_vector2];
end

%6.5计算reuters的平均值
% 初始化新向量 prethetare
prethetare = [];
prevertare = [];

% 循环计算每个向量的平均值
for i = 1:6
    % 获取当前向量
    current_vector = discomre.(['disrea', num2str(i)]);  % 假设向量保存在结构体的 vector 字段中
    current_vector2 = discomre.(['disreb', num2str(i)]);  
    % 计算当前向量的平均值
    mean_vector = mean(current_vector);
    mean_vector2 = mean(current_vector2);
    % 将平均值追加到向量 prethetabb 中
    prethetare = [prethetare, mean_vector];
    prevertare = [prevertare, mean_vector2];
end

%6.6计算guardian的平均值
% 初始化新向量 prethetagu
prethetagu = [];
prevertagu = [];

% 循环计算每个向量的平均值
for i = 1:6
    % 获取当前向量
    current_vector = discomgu.(['disgua', num2str(i)]);  % 假设向量保存在结构体的 vector 字段中
    current_vector2 = discomgu.(['disgub', num2str(i)]);  
    % 计算当前向量的平均值
    mean_vector = mean(current_vector);
    mean_vector2 = mean(current_vector2);
    % 将平均值追加到向量 prethetabb 中
    prethetagu = [prethetagu, mean_vector];
    prevertagu = [prevertagu, mean_vector2];
end


%      7 使用维度转化公式对数据进行处理，将bbcpretheta1-6和bbcpreverta1-6转化为bbctheta1-6 和 bbcverta1-6
%现在，需要设定维度转的话参数 k_dimention=0.95到1.00对数据进行处理。
kdim=0.97;  

%！！！！重要！！参数设定！！
thetabb= kdim * prethetabb;
vertabb= kdim * prevertabb;
thetare= kdim * prethetare;
vertare= kdim * prevertare;
thetagu= kdim * prethetagu;
vertagu= kdim * prevertagu;


%      8 对theta1-6 和 verta1-6 进行处理得到热核函数的两个重要参数 sigma1-6 和 eta1-6
sigmabb= repelem(thetabb,1,11);
etabb= repelem( (vertabb / (1.7338*sqrt(2))),1,11);
sigmare= repelem(thetare,1,11);
etare= repelem( (vertare / (1.7338*sqrt(2))),1,11);
sigmagu= repelem(thetagu,1,11);
etagu= repelem( (vertagu / (1.7338*sqrt(2))),1,11);


%      9 对于bag1-6,计算bag中所有样本到cMbbcA1-6的距离再用维度转换公式得到子类所有bbcZ_i

%    9.1 将discombb这个结构体数组中名字带有disbbtotal1到disbbtotal6这几个向量进行拼接用作Z_i的初始化
% 初始化一个空向量，用于存储拼接后的结果
convecbb = [];
% 循环遍历 discombb 结构体数组，拼接符合条件的向量
for i = 1:6
    field_name = ['disbbtotal', num2str(i)];
    current_vector = discombb.(field_name);
    convecbb = [convecbb; current_vector];
end
convecbb =convecbb' ;
%    9.2 对reuters数据进行处理
convecre = [];
% 循环遍历 discombb 结构体数组，拼接符合条件的向量
for i = 1:6
    field_name = ['disretotal', num2str(i)];
    current_vector = discomre.(field_name);
    convecre = [convecre; current_vector];
end
convecre =convecre' ;
%    9.3 对guardian数据进行再次处理
convecgu = [];
% 循环遍历 discombb 结构体数组，拼接符合条件的向量
for i = 1:6
    field_name = ['disbbtotal', num2str(i)];
    current_vector = discombb.(field_name);
    convecgu = [convecgu; current_vector];
end
convecgu =convecgu' ;


%      10 计算
%      C_ii=exp(-(bbcZ_i)^2/(2bbctheta1-6^2))+exp(-(bbcZ_i)^2/(2bbcverta1-6^2))?????写错了？？？
%10.1 计算bbc上的C矩阵向量

%设置权重参数gamma，默认eta的系数小一些
gamma=0.25;

%cbb = exp(-(convecbb.^2)./(2*sigmabb.^2)) + exp(-(convecbb.^2)./(2*etabb.^2));
%cre = exp(-(convecre.^2)./(2*sigmare.^2)) + exp(-(convecre.^2)./(2*etare.^2));
%cgu = exp(-(convecgu.^2)./(2*sigmagu.^2)) + exp(-(convecgu.^2)./(2*etagu.^2));

cbb2= gamma*(exp(-(convecbb.^2)./(2*etabb.^2)))./etabb;
cbb1= (1-gamma)*(exp(-(convecbb.^2)./(2*sigmabb.^2)))./sigmabb;
cbb2(isnan(cbb2)) = 0;
cbb=cbb1+cbb2;

cre2= gamma*(exp(-(convecre.^2)./(2*etare.^2)))./etare;
cre1= (1-gamma)*(exp(-(convecre.^2)./(2*sigmare.^2)))./sigmare;
cre2(isnan(cre2)) = 0;
cre=cre1+cre2;

cgu2= gamma*(exp(-(convecgu.^2)./(2*etagu.^2)))./etagu;
cgu1= (1-gamma)*(exp(-(convecgu.^2)./(2*sigmagu.^2)))./sigmagu;
cgu2(isnan(cgu2)) = 0;
cgu=cgu1+cgu2;


Cbb=diag(100.*cbb);
Cre=diag(100.*cre);
Cgu=diag(100.*cgu);

%527行

% 假设 indicatorguardian 是一个包含 6 个结构体的结构体数组


% 初始化一个空向量，用于存储连接后的结果
outliers_x = [];

% 循环访问每个结构体的字段，并将其连接到结果向量中
for i = 1:6
    field_name = sprintf('ing%d', i);  % 构建字段名称，如 'ing1' 到 'ing6'
    current_field = indicatorguardian.(field_name);  % 获取当前结构体的字段值
    outliers_x = [outliers_x, current_field];  % 将字段值连接到结果向量中
end

% 显示结果


