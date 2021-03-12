function [sortedEigVecB,singularValueMatrix,sortedEigVecA]=mysvd(matrixToCompute)
    %% 形成左右矩阵。其中A对应右侧V，B对应左侧U
    matrixA = matrixToCompute.'*matrixToCompute;
    matrixB = matrixToCompute*matrixToCompute.';
    %% 记待求矩阵为A，求解A'A和AA'
    [eigVecA,eigValueA]=eig(matrixA);
    [eigVecB,eigValueB]=eig(matrixB);
    %% 对特征值和特征向量做排序
    % 排序算法来自 https://ww2.mathworks.cn/help/matlab/ref/eig.html
    [eigValueElementA,sortedIndexA] = sort(diag(eigValueA),'descend');
    sortedEigValA=eigValueA(sortedIndexA,sortedIndexA);
    sortedEigVecA = eigVecA(:,sortedIndexA);
    
    %% 对 matrixB 的特征值和特征向量做同样的排序操作
    [eigValueElementB,sortedIndexB] = sort(diag(eigValueB),'descend');
    sortedEigValB=eigValueB(sortedIndexB,sortedIndexB);
    sortedEigVecB = eigVecB(:,sortedIndexB);
    %% 截取有效部分
    if length(sortedEigValB)>length(sortedEigValA)
        singularValueMatrix=sqrt(sortedEigValB(1:size(matrixToCompute,1),1:size(matrixToCompute,2)));
    else
        singularValueMatrix=sqrt(sortedEigValA(1:size(matrixToCompute,1),1:size(matrixToCompute,2)));
    end
    %% 检查左右U,V是否能够满足中间的奇异值分解主对角线元素非负。不满足则修改左侧使之满足
    for i= 1: length(find(eigValueA>1e-6))
        if max(abs((1/sqrt(sortedEigValA(i,i))*matrixToCompute*sortedEigVecA(:,i) + sortedEigVecB(:,i))))<1e-6
            sortedEigVecB(:,i)=-sortedEigVecB(:,i);
        end    
    end
    %% 验证SVD是否成功
    sortedEigVecB*singularValueMatrix*sortedEigVecA.';
end