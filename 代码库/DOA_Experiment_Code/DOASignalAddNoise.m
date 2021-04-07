function noiseMatrix = DOASignalAddNoise(Ps,SNR,arrayElementNumber,sampleNumber)
    % ��ģ�����AWGN�ŵ�ģ�ͣ�������Դ�����źŵĹ��ʺ�ʵ����ָ��������ȣ�
    % ����һ����˹����������
    
    % ���������������壺
    % Ps���źŵķ��书�ʣ����ź�����ģ�鷵�ء�
    % SNR������ȣ���ʵ����ָ����
    % arrayElementNumber����Ԫ��������DOA����ʵ���ȫ�ֱ�����
    % sampleNumber��ʱ������ĵ�������DOA����ʵ���ȫ�ֱ�����
    
    % �����������Ϊ��������
    noiseMatrix = (Ps/10^SNR)... % ��������ȼ���õ�����ǿ��
    *(randn(arrayElementNumber,sampleNumber)+randn(arrayElementNumber,sampleNumber)*1i);
    % ���ɸ�˹��������    
end