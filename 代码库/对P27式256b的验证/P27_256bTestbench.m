clear;clc;
% 可配置项，测试的总次数
testCount=10000;
% 可配置项，测试的误差容限
tolerableErr = 1e-11;

% 生成矩阵的随机维度
randomDim=fix(200*rand(1,testCount))+2;

% 接收矩阵初始化
sizeInformation = [];

for index=1:testCount
    % 生成随机维度的 R 矩阵
    [signalDim,R]=RMatGen(randomDim(index));
    % 求 R 矩阵的特征值和特征向量
    [U,Sigma]=eig(R);
    % 对求得的特征值矩阵降序排序
    [~,sortedIndex] = sort(diag(Sigma),'descend');
    % 以对应的顺序对特征向量排序
    sortedSigma = Sigma(sortedIndex,sortedIndex);
    sortedU = U(:,sortedIndex);
    
    % 提取 Us 和 Un
    Usignal = sortedU(:,1:signalDim);
    Unoise = sortedU(:,signalDim+1:end);
    
    % 验证 2.5.6a
    answerOf256a=Usignal*Usignal'+Unoise*Unoise';
    % 如果测试用例的个数不超过 10 个，
    % 记录 2.5.6a 的答案，即Usignal*Usignal'+Unoise*Unoise'矩阵，的尺寸与当前 R 矩阵的尺寸
    % 通过观察此列表，可以得出 2.5.6a 的答案是和R同阶的矩阵的结论。
    % 限制测试用例数量的原因是当测试数量很大的时候，变长数组需要的计算机硬件性能较高，生成这个数组需要较长的时间。
    if(testCount<=10)
        sizeInformation = [sizeInformation;randomDim(index) size(answerOf256a)];
    end
    % 将 2.5.6a 的答案与同阶的单位矩阵作差，求其误差；
    % 然后与误差容限比较，如果超出误差容限，errorInformation 矩阵就会获得元素从而非空。
    % 若 errorInformation 为空矩阵，则可以验证 2.5.6a 成立。
    errorInformation = find((answerOf256a-eye(randomDim(index)))>tolerableErr);
    
    % 验证 2.5.6b
    UsUsH = Usignal*Usignal';
    UnUnH = Unoise*Unoise';
    % UnUnH 的非零元都集中在右下角不便操作，将其元素降序排序，使其元素移到左上角
    [~,sortedUnUnHIndex] = sort(diag(UnUnH),'descend');
    sortedUnUnH = UnUnH(sortedUnUnHIndex,sortedUnUnHIndex);
    % 将 Us*Us' 与 Un*Un' 取其非零部分，将取出的非零部分与和它同秩的单位矩阵作差，求其误差；
    % 然后与误差容限比较，如果超出误差容限，errorInformationUsUsH 和 errorInformationUnUnH 矩阵就会获得元素从而非空。
    % 若 errorInformationUsUsH 和 errorInformationUnUnH 矩阵均为空矩阵，则可以验证 2.5.6b 成立。
    errorInformationUsUsH = find((UsUsH(1:rank(UsUsH),1:rank(UsUsH))-eye(rank(UsUsH)))>tolerableErr);
    errorInformationUnUnH = find((UnUnH(1:rank(UnUnH),1:rank(UnUnH))-eye(rank(UnUnH)))>tolerableErr);
end
