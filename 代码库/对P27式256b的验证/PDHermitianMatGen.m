function PDHermitianMatrix=PDHermitianMatGen(dim)
    % 该函数用于生成正定Hermitian矩阵
    % 由于27页的式2.5.6b所用的信号子空间矩阵Us和噪声子空间矩阵Un，
    % 均来自正定Hermitian矩阵R的酉矩阵对角化，
    % 这里用这个函数来生成正定Hermitian矩阵作为矩阵R。
  
    % 正定Hermitian矩阵生成算法
    PDHermitianMatrix=rand(dim)+i*rand(dim);
    PDHermitianMatrix=PDHermitianMatrix'*PDHermitianMatrix+0.01*eye(dim);
end
