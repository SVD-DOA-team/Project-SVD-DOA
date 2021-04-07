% DOA仿真脚本

% 全局变量设置部分
% 设置信源的个数
signalNumber = 5;
% 设置阵元的个数，需要改变参数时，把默认值4改成其他的数值即可。
% 这里使用了max()，保证阵元个数至少比信源个数多1。
arrayElementNumber = max(20, signalNumber+1);
% 设置时域采样的点数。
sampleNumber = 120;
% 真实的各波达角组成的向量（角度制，非弧度制）
% 既是DOA估计的参考答案，也是阵列流形生成算法所需的参数。
realDOA = [-80 -40 15 60 75];

% 仿真实验算法
% 信源发出信号
% 信源生成算法的可动参数：
% 载频f
f = 13.56*10^6;
% 载频幅度
amplitude = 10;
[lambda,Ps,signalsOutput]=DOASignalGenerating(signalNumber,sampleNumber,f,amplitude);
% 信道的加性高斯噪声
% 信道噪声生成算法的可动参数：
% 信噪比SNR
SNR = 1.5;
noiseMatrix = DOASignalAddNoise(Ps,SNR,arrayElementNumber,sampleNumber);
% 阵列接收信号和处理
% 信号接收和处理部分的可动参数：
% 阵元间隔d
d = lambda/2;
DOASignalProceeding(signalNumber,d,lambda,arrayElementNumber,signalsOutput,noiseMatrix,realDOA,sampleNumber);
