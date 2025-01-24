
% 生成1到8的整数序列
sequence = 1:8;

% 随机排列整数序列
random_sequence = randperm(length(sequence));

% 从随机序列中选择前5个数
selected_numbers = sequence(random_sequence(1:5));

disp('随机选择的5个数：');
disp(selected_numbers);