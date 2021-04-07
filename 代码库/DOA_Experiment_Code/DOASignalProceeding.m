function DOASignalProceeding...
    (signalNumber,...
    arrayEleDistance,...
    lambda,...
    arrayElementNumber,...
    signalsOutput,...
    noiseMatrix,...
    realDOA,...
    sampleNumber)

    % 本模块代表信号接收端，接收发送端发送的单音信号，估计其波达角。
    
    % 其各输入参数的意义：
    % signalNumber是信号个数。
    % arrayEleDistance是阵元间距，这个参数的取值与空域采样的奈奎斯特定理有关。
    % 是一个由实验者自行设定的参数。
    % lambda是载波波长。
    % arrayElementNumber是阵元个数，是DOA仿真实验的全局变量。
    % signalsOutput是信号源发送的信号。
    % noiseMatrix是信道产生的加性高斯噪声。
    % 这两个参数用于与阵列流形合成接收到的信号。
    % realDOA是真实的波达角，是DOA仿真实验的全局变量。
    % 在这里我们需要知道真实的波达角来生成阵列流形。
    % sampleNumber是时域采样的点数，是DOA仿真实验的全局变量。
 
    % 生成阵列流形
    a0=-1i*2*pi*arrayEleDistance/lambda.*(0:arrayElementNumber-1);
    A0=exp(a0.'*sin(dec2rad(realDOA)));

    %阵列输出信号采样
    x=A0*signalsOutput+noiseMatrix;

    %求出采样R矩阵
    R=(x*x')./sampleNumber;

    % 对R进行特征分解，得到U，Un，Us
    % 这里特意使用更高效的svd()函数，
    % 数学上可以证明这里使用svd()和eig()等效
    [U,~]=svd(R);
    Un=U(:,signalNumber+1:arrayElementNumber);

    % 设扫描精度为L
    L=500;
	% 让theta在-90度到90度之间变化
    theta=-90:180/L:90-180/L;
    % 依据公式计算P(\theta)
    a=exp(-1i.*2.*pi.*arrayEleDistance./lambda*(0:arrayElementNumber-1)'*sin(dec2rad(theta)));
    P=diag(1./(a'*(Un*Un')*a));

    % 画图
    plot(theta,10*log10(abs(P/max(P))));grid;
    axis([-90,90,-20,0]);
    xlabel('\theta/degree');ylabel('P(\theta)/dB');
    title('P(\theta)');
end