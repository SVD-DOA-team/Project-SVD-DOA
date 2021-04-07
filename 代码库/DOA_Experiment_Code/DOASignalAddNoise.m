function noiseMatrix = DOASignalAddNoise(Ps,SNR,arrayElementNumber,sampleNumber)
    % 本模块仿真AWGN信道模型，根据信源发射信号的功率和实验者指定的信噪比，
    % 产生一个高斯白噪声矩阵。
    
    % 其各输入参数的意义：
    % Ps是信号的发射功率，由信号生成模块返回。
    % SNR是信噪比，由实验者指定。
    % arrayElementNumber是阵元个数，是DOA仿真实验的全局变量。
    % sampleNumber是时域采样的点数，是DOA仿真实验的全局变量。
    
    % 其输出参数即为噪声矩阵。
    noiseMatrix = (Ps/10^SNR)... % 根据信噪比计算得到噪声强度
    *(randn(arrayElementNumber,sampleNumber)+randn(arrayElementNumber,sampleNumber)*1i);
    % 生成高斯噪声矩阵    
end