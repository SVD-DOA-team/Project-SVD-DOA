function [sortedEigVecU,singularValueMatrix,sortedEigVecV]=homemade_sorted_svd(matrixToCompute)
    % 将转置.'改为共轭转置'，使其支持复数矩阵输入
    baseMatrixV = matrixToCompute'*matrixToCompute;
    baseMatrixU = matrixToCompute*matrixToCompute';
    % 记待求矩阵为A，求解A'A和AA'
    [eigVecV,eigValueV]=eig(baseMatrixV);
    [eigVecU,eigValueU]=eig(baseMatrixU);
    % 对特征值和特征向量做排序
    % 排序算法来自 https://ww2.mathworks.cn/help/matlab/ref/eig.html
    % 排序的目的是方便之后裁剪奇异值矩阵
    [~,sortedIndexV] = sort(diag(eigValueV),'descend');
    % 抛掉不要的特征值向量，优化程序性能
    sortedEigValV = eigValueV(sortedIndexV,sortedIndexV);
    sortedEigVecV = eigVecV(:,sortedIndexV);
    
    [~,sortedIndexU] = sort(diag(eigValueU),'descend');
    % 抛掉不要的特征值向量，优化程序性能
    sortedEigValU=eigValueU(sortedIndexU,sortedIndexU);
    sortedEigVecU = eigVecU(:,sortedIndexU);
    vectorOfEigVal = sqrt(diag(sortedEigValU));
    
    % 纠正左奇异矩阵的各列向量的符号
    for index=1:min(length(sortedEigVecU),length(sortedEigVecV))
    % index 的取值范围取决于左右奇异矩阵中尺寸较小的一个的尺寸
        if abs(matrixToCompute*sortedEigVecV(:,index)+vectorOfEigVal(index)*sortedEigVecU(:,index))<1e-6
        % 针对复数的情况，改为判断复数的模
        % 左右奇异矩阵的列向量必须满足 A*vi=sigma*ui 这一等量关系
        % 如果相差一个负号，必有 A*vi=-sigma*ui
        % 为了浮点数检验方便，将上式移项，改为检查 A*vi+sigma*ui<误差上限
            sortedEigVecU(:,index)=-sortedEigVecU(:,index);
            % 若满足“相差一个负号”的条件，就将这一列列向量反号
        end
    end
    % 将特征值矩阵裁剪为奇异值矩阵的算法
    % 选取两个特征值矩阵中尺寸较大的一个，取其前m行、前n列（若输入矩阵为m*n)
    if length(sortedEigValU)>length(sortedEigValV)
        singularValueMatrix=sqrt(sortedEigValU(1:size(matrixToCompute,1),1:size(matrixToCompute,2)));
    else
        singularValueMatrix=sqrt(sortedEigValV(1:size(matrixToCompute,1),1:size(matrixToCompute,2)));
    end
end
