%%迭代过程设计《奥米加档案-极昼》
%第一板块中存在一些特殊参数现在需要再进行确认一下
%buildingC中的几个参数
%1.维度转换参数kdim,和视角无关
kdim=0.95;
%2.gamma 多热核函数的均衡参数
gamma=0.25;
%3.alpha0 损失函数和图正则化项的均衡参数
alpha0=0.1;
%4.beta0 沟通项和其他项的均衡参数
beta0=1;
%除此之外惩罚图和超图之间的平衡参数应该如何设定？？？
%图Laplace矩阵的混合比率按照2比1进行混合
%在横剖其
%平衡系数沟通
beta1=1.4113;
beta2=1.3238;
beta3=1.7115;
%由于范数平方的原因，所以需要将范数进行平方
preBeta1=1/beta1^2;
preBeta2=1/beta2^2;
preBeta3=1/beta3^2;
sumpreBeta=preBeta1+preBeta2+preBeta3;

Beta1=preBeta1/sumpreBeta;
Beta2=preBeta2/sumpreBeta;
Beta3=preBeta3/sumpreBeta;


%



%全1矩阵

%准备图Laplace矩阵我们只用L来进行表示，其中LH1和Lh分别代表两种不同性质的图
LH1=newlaplacebb;
LH2=newlaplacere;
LH3=newlaplacegu;

Lh1=newlaplace2bb;
Lh2=newlaplace2re;
Lh3=newlaplace2gu;

L1=2.*LH1+Lh1;
L2=2.*LH2+Lh2;
L3=2.*LH3+Lh3;
%拆分为positive和minus部分
%L1拆分
diag_elements1 = diag(L1);
L1positive = diag(diag_elements1);
L1minus = -1*(L1- L1positive);
%L2拆分
diag_elements2 = diag(L2);
L2positive = diag(diag_elements2);
L2minus = -1*(L2- L2positive);
%L3拆分
diag_elements3 = diag(L3);
L3positive = diag(diag_elements3);
L3minus = -1*(L3- L3positive);



%%以下是一个C的预备
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

C1_old=Cbb;
C2_old=Cre;
C3_old=Cgu;

X1=m_scaled_bbc'+0.05;%原始数据矩阵记为X，后第一位数字代表模态
X2=m_scaled_reuters'+0.05;
X3=m_scaled_guardian'+0.05;

%对于迭代过程中会发生变化的量拆分成两代 
B1_old=B1+0.05;
B2_old=B2+0.05;
B3_old=B3+0.05;

Q1=diag(sum(B1_old));
%计算B2对应的Q2
Q2=diag(sum(B2_old));
%计算B3对应的Q3
Q3=diag(sum(B3_old));
Q1_old=Q1;
Q2_old=Q2;
Q3_old=Q3;

F1=PreF1;%F的初始化
F2=PreF2;
F3=PreF3;
F1_old=F1;
F2_old=F2;
F3_old=F3;

Fx = Beta1*PreF1*Q1+Beta2*PreF2*Q2+Beta3*PreF3*Q3; 
Fx_old=Fx;

M11= ones(3560, 6);
M12= ones(3068, 6);
M13= ones(3631, 6);
%现在需要额外计算R矩阵和S矩阵



for i=1:1
i=i+1
%A1   更新策略以及更新过程
R1positive=F1_old'*F1_old*Q1_old;
R1minus=F1_old'*Fx_old;
S1positive=F1_old'*L1positive*F1_old*Q1_old;
S1minus=F1_old'*L1minus*F1_old*Q1_old;
%更新B并且将B相关的q也进行更新,并将更新后的矩阵进行重新赋值和储存

%R是沟通，S是图
%检查数量级
%btest1=X1*C1_old*F1_old 
%btest2=Beta1*M11*R1minus
%btest3=alpha0*M11*S1minus

B1up= X1*C1_old*F1_old + Beta1*M11*R1minus+alpha0*M11*S1minus;
B1down=B1_old*F1_old'*C1_old*F1_old+Beta1*M11*R1positive+alpha0*M11*S1positive;
B1_new=B1up./B1down;

B1_old=B1_new;
Q1_old=diag(sum(B1_new));
%检查f数量级
%ftest1=C1_old*X1'*B1_old; 
%ftest2= Beta1*F1_old * Q1_old;
%ftest3=2*alpha0*L1minus*F1_old*Q1_old*Q1_old;
%fdtest1=C1_old*F1_old*B1_old'*B1_old;
F1up= C1_old*X1'*B1_old + Beta1*Fx_old * Q1_old + alpha0*L1minus*F1_old*Q1_old*Q1_old;
F1down= C1_old*F1_old*B1_old'*B1_old + Beta1*F1*Q1_old*Q1_old + alpha0*L1positive*F1_old*Q1_old*Q1_old;
F1_new=F1up./F1down;

F1_old=F1_new

%更新F之后，我们需要计算%X-BF%进而计算C1我们值得注意的是
Di_new = X1-B1_old*F1_old';
vecnormX_new1 = vecnorm(Di_new);
c12=gamma * (exp(-(vecnormX_new1 .^2)./(2*etabb.^2))) ./ etabb;
c11=(1-gamma)* (exp(-(vecnormX_new1 .^2)./(2*sigmabb.^2))) ./ sigmabb;
c12(isnan(c12)) = 0;
c1=c11+c12;
C1_old=diag(100.*c1);%这个系数可能徐涛调整


%以上C1更新完毕
R2positive=F2_old'*F2_old*Q2_old;
R2minus=F2_old'*Fx_old;
S2positive=F2_old'*L2positive*F2_old*Q2_old;
S2minus=F2_old'*L2minus*F2_old*Q2_old;
%更新B并且将B相关的q也进行更新,并将更新后的矩阵进行重新赋值和储存
B2up= X2*C2_old*F2_old + Beta2*M12*R2minus + alpha0*M12*S1minus;
B2down=B2_old*F2_old'*C2_old*F2_old + Beta2*M12*R2positive + alpha0*M12*S2positive;
B2_new=B2up./B2down;
B2_old=B2_new;
Q2_old=diag(sum(B2_new));
%更新F
F2up=C2_old*X2'* B2_old  +  Beta2 * Fx_old * Q2_old  +  alpha0*L1minus*F2_old*Q2_old*Q2_old;
F2down=C2_old*F2_old*B2_old'*B2_old + Beta2*F2*Q2_old*Q2_old + alpha0*L2positive*F2_old*Q2_old*Q2_old;
F2_new=F2up./F2down;
F2_old=F2_new;
%更新F之后，我们需要计算%X-BF%进而计算C1我们值得注意的是
Di_new = X2-B2_old*F2_old';
vecnormX_new2 = vecnorm(Di_new);
c22=gamma*(exp(-(vecnormX_new2 .^2)./(2*etabb.^2)))./etabb;
c21=(1-gamma)*(exp(-(vecnormX_new2 .^2)./(2*sigmabb.^2)))./sigmabb;
c22(isnan(c22)) = 0;
c2=c21+c22;
C2_old=diag(100.*c1);
%以上C2更新完毕

R3positive=F3_old'*F3_old*Q3_old;
R3minus=F3_old'*Fx_old;
S3positive=F3_old'*L3positive*F3_old*Q3_old;
S3minus=F3_old'*L3minus*F3_old*Q3_old;
%更新B并且将B相关的q也进行更新,并将更新后的矩阵进行重新赋值和储存
B3up= X3*C3_old*F3_old + Beta3*M13*R3minus + alpha0*M13*S3minus;
B3down=B3_old*F3_old'*C3_old*F3_old + Beta3*M13*R3positive + alpha0*M13*S3positive;
B3_new=B3up./B3down;
B3_old=B3_new;
Q3_old=diag(sum(B3_new));
%更新F
F3up=C3_old*X3'*B3_old + Beta3*Fx_old*Q3_old + alpha0*L3minus*F3_old*Q3_old*Q3_old;
F3down=C3_old*F3_old*B3_old'*B3_old + Beta3*F3*Q3_old*Q3_old + alpha0*L3positive*F3_old*Q3_old*Q3_old;
F3_new=F3up./F3down;
F3_old=F3_new;
%更新F之后，我们需要计算%X-BF%进而计算C1我们值得注意的是
Di_new = X3-B3_old*F3_old';
vecnormX_new3 = vecnorm(Di_new);
c32=gamma*(exp(-(vecnormX_new3 .^2)./(2*etabb.^2)))./etabb;
c31=(1-gamma)*(exp(-(vecnormX_new3 .^2)./(2*sigmabb.^2)))./sigmabb;
c32(isnan(c32)) = 0;
c3=c31+c32;
C3_old=diag(100.*c3);

%以上C3更新完毕
Fx_old = Beta1*F1_old*Q1_old + Beta2*F2_old*Q2_old  + Beta3*F3_old*Q3_old 
end



alpha1=1;
alpha0=1;N11=ones(6,6);
%备选更新方案，将模态分解为多个独立部分再组合,迈进，迈向前方
for i=1:50
%C1_old=eye(66);

B1up= X1*C1_old*F1_old + 1000*B1_old;
B1down=B1_old*F1_old'*C1_old*F1_old + 1000*B1_old*N11';

B1_old=B1up./B1down;

F1up=C1_old*X1'*B1_old + Beta1*F1_old +0.001* L1minus*F1_old;
F1down=C1_old*F1_old*B1_old'*B1_old + Beta1*F1_old + 0.001*L1positive*F1_old;

F1_old=F1up./F1down;
Di_new = X1-B1_old*F1_old';
vecnormX_new1 = vecnorm(Di_new);
c12=gamma * (exp(-(vecnormX_new1 .^2)./(2*etabb.^2))) ./ etabb;
c11=(1-gamma)* (exp(-(vecnormX_new1 .^2)./(2*sigmabb.^2))) ./ sigmabb;
c12(isnan(c12)) = 0;
c1=c11+c12;
C1_old=diag(100.*c1);
end






%2计算F1的更新；更新后记载并储存

%3计算B2的更新；计算Q2（B2）的更新；更新后记载并储存

%4计算F2的更新；更新后记载并储存

%5计算B3的更新；计算Q3（B3）的更新；更新后记载并储存

%6计算F3的更新；更新后记载并储存

%7更新C1矩阵；更新后记载并储存

%8更新C2矩阵；更新后记载并储存

%9更新C3矩阵；更新后记载并储存

%10更新F*，或者称之为Fx

