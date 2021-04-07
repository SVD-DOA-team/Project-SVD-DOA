function [lambda,Ps,signalsOutput]=DOASignalGenerating...
    (signalNumber,...
    sampleNumber,...
    freqCarrier,...
    ampSignal)
    % ��ģ������źŷ��Ͷˣ�����һϵ�е����źš�
    
    % ���������������壺
    % ������������������Դ�������ز�Ƶ�ʡ��źŷ��ȣ��������źš�
    % ����������ָ�����ɵĸ������źŽ���ʱ������ĵ�����
    % �������źž��������������������������һ��ȫ�ֱ�����
    
    % ���������������壺
    % lambda��ʾ�ز��Ĳ������������ز�Ƶ�ʡ�
    % Ps���źŵķ��书�ʣ������ŵ�ģ��������������������������Ȳ���������������
    % signalsOutput��K����Դ���͵�K�������ź���ɵľ���
    
    % ���㲢�������
    % ����
    c=3*10^8;
    % ����õ�����
    lambda=c/freqCarrier;
    
    % ���㲢����źŷ��书��
    % ����Ĺ��������ҵ����źŵľ������ʡ�
    % ����źŲ������ҵ����źţ�����Ҫ�޸Ĺ��ʵĹ�ʽ��
    Ps = ampSignal^2/2;
    
    % ���㲢����źž���
    % ���ɷ����źŵĽ�Ƶ��
    w = (pi/2/signalNumber).*(1:signalNumber);
    % ���ɲ���ʱ������
    t = 0:sampleNumber-1;
    % �����źž���
    signalsOutput = zeros(signalNumber,sampleNumber);
    for index = 1:signalNumber
    signalsOutput(index,:)=ampSignal*cos(w(index)*t);
    end
end