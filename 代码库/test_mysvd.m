clc;clear;close all;
%% 可配置项，修改此处来决定可容忍的误差精度
tolerableErr = 1e-8;
%% 初始化保存运算误差的两个变量
errS=[];
err=[];
%% 循环检验算法的有效性
for i=1:10000
   a=fix(200*rand)+2;%随机行数
   b=fix(200*rand)+2;%随机列数
   mat=rand(a,b)+i*rand(a,b);%随机矩阵
   [stdU,stdS,stdV]=svd(mat);%标准参考
   [myU,myS,myV]=homemade_sorted_svd(mat);%自己的运算结果
   errS=[errS,max(abs(myS-stdS))];%奇异值的差
   err=[err,max(max(abs(myU*myS*myV'-mat)))];%还原之后相比于原矩阵的差
end
%% 输出最大值
maxErrS=max(errS)
maxErr=max(err)

%% 输出算法出错的频率
freqErrS=length(find(errS>tolerableErr))/length(errS)
freqErr=length(find(err>tolerableErr))/length(err)
