% DOA����ű�

% ȫ�ֱ������ò���
% ������Դ�ĸ���
signalNumber = 5;
% ������Ԫ�ĸ�������Ҫ�ı����ʱ����Ĭ��ֵ4�ĳ���������ֵ���ɡ�
% ����ʹ����max()����֤��Ԫ�������ٱ���Դ������1��
arrayElementNumber = max(20, signalNumber+1);
% ����ʱ������ĵ�����
sampleNumber = 120;
% ��ʵ�ĸ��������ɵ��������Ƕ��ƣ��ǻ����ƣ�
% ����DOA���ƵĲο��𰸣�Ҳ���������������㷨����Ĳ�����
realDOA = [-80 -40 15 60 75];

% ����ʵ���㷨
% ��Դ�����ź�
% ��Դ�����㷨�Ŀɶ�������
% ��Ƶf
f = 13.56*10^6;
% ��Ƶ����
amplitude = 10;
[lambda,Ps,signalsOutput]=DOASignalGenerating(signalNumber,sampleNumber,f,amplitude);
% �ŵ��ļ��Ը�˹����
% �ŵ����������㷨�Ŀɶ�������
% �����SNR
SNR = 1.5;
noiseMatrix = DOASignalAddNoise(Ps,SNR,arrayElementNumber,sampleNumber);
% ���н����źźʹ���
% �źŽ��պʹ����ֵĿɶ�������
% ��Ԫ���d
d = lambda/2;
DOASignalProceeding(signalNumber,d,lambda,arrayElementNumber,signalsOutput,noiseMatrix,realDOA,sampleNumber);
