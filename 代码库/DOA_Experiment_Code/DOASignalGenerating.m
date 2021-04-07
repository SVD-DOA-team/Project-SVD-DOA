function [lambda,Ps,signalsOutput]=DOASignalGenerating...
    (signalNumber,...
    sampleNumber,...
    freqCarrier,...
    ampSignal)
    % 本模块代表信号发送端，生成一系列单音信号。
    
    % 其各输入参数的意义：
    % 接收三个参数——信源个数、载波频率、信号幅度，来生成信号。
    % 采样个数是指对生成的各正弦信号进行时域采样的点数，
    % 将决定信号矩阵、噪声矩阵的列数，本质上是一个全局变量。
    
    % 其各输出参数的意义：
    % lambda表示载波的波长，反比于载波频率。
    % Ps是信号的发射功率，仿真信道模块接受这个参数，并依据其信噪比参数来生成噪声。
    % signalsOutput是K个信源发送的K个基带信号组成的矩阵。
    
    % 计算并输出波长
    % 光速
    c=3*10^8;
    % 计算得到波长
    lambda=c/freqCarrier;
    
    % 计算并输出信号发射功率
    % 这里的功率是正弦单音信号的均方功率。
    % 如果信号不是正弦单音信号，可能要修改功率的公式。
    Ps = ampSignal^2/2;
    
    % 计算并输出信号矩阵
    % 生成发射信号的角频率
    w = (pi/2/signalNumber).*(1:signalNumber);
    % 生成采样时刻向量
    t = 0:sampleNumber-1;
    % 生成信号矩阵
    signalsOutput = zeros(signalNumber,sampleNumber);
    for index = 1:signalNumber
    signalsOutput(index,:)=ampSignal*cos(w(index)*t);
    end
end