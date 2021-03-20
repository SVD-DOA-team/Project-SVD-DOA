clear;clc;
% 可配置项，测试的总次数
testCount=10000;
PDflags=[];
HMflags=[];
randomDim=fix(200*rand(1,testCount))+2;
for index=1:testCount
    mat=PDHermitianMatGen(randomDim(index));
    % 判断是否为Hermitian矩阵
    isHM=ishermitian(mat);
    % 判断是否为正定矩阵
    try chol(mat);
        isPD=1;
    catch ME
        isPD=0;
    end
    HMflags=[HMflags isHM];
    PDflags=[PDflags isPD];
end
indexOfNotPD=find(~PDflags)
indexOfNotHM=find(~HMflags)
  
% 这一测试脚本会调用PDHermitianMatGen函数，生成一定数量的2-202阶待检测方阵。
% 对生成的矩阵mat判断是否Hermitian和是否正定，并记录判断结果。
% 最后用find函数找出不满足Hermitian或正定的矩阵所在的索引。
% 如果indexOfNotPD 和 indexOfNotHM 均为空矩阵，
% 则 PDHermitianMatGen 函数没有出错。
