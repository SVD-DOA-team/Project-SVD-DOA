function [signalDim,RMatrix]=RMatGen(dim)
    % 该函数用于生成满足公式需要的R矩阵
    % 27页的式2.5.6b所用的信号子空间矩阵Us和噪声子空间矩阵Un，
    % 均来自R矩阵的酉矩阵对角化，
    % R矩阵不仅正定、Hermitian，而且其最小的特征值重复若干次（即噪声分量）。
    % 这使得R矩阵的特征值可以分为两部分，构成了Us和Un的划分基础。
    % 因此我们用更复杂的方法生成R矩阵。
  
    % 决定信号分量的维数
    signalDim=ceil((dim-1)*rand);
    
    % 决定噪声功率sigma^2
    noiseCoefficent=rand;
    
    % 构造R矩阵
    % 构造信号部分，一个半正定的Hermitian矩阵
    RMatrix=rand(signalDim)+1i*rand(signalDim);
    RMatrix=RMatrix'*RMatrix;
    % 拓展R矩阵的维度，为噪声分量留出空间
    RMatrix(dim,dim)=0;
    % 加入噪声扰动，同时使得R正定
    RMatrix=+noiseCoefficent*eye(dim);
end
