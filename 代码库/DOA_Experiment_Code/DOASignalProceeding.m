function DOASignalProceeding...
    (signalNumber,...
    arrayEleDistance,...
    lambda,...
    arrayElementNumber,...
    signalsOutput,...
    noiseMatrix,...
    realDOA,...
    sampleNumber)

    % ��ģ������źŽ��նˣ����շ��Ͷ˷��͵ĵ����źţ������䲨��ǡ�
    
    % ���������������壺
    % signalNumber���źŸ�����
    % arrayEleDistance����Ԫ��࣬���������ȡֵ�����������ο�˹�ض����йء�
    % ��һ����ʵ���������趨�Ĳ�����
    % lambda���ز�������
    % arrayElementNumber����Ԫ��������DOA����ʵ���ȫ�ֱ�����
    % signalsOutput���ź�Դ���͵��źš�
    % noiseMatrix���ŵ������ļ��Ը�˹������
    % �����������������������κϳɽ��յ����źš�
    % realDOA����ʵ�Ĳ���ǣ���DOA����ʵ���ȫ�ֱ�����
    % ������������Ҫ֪����ʵ�Ĳ�����������������Ρ�
    % sampleNumber��ʱ������ĵ�������DOA����ʵ���ȫ�ֱ�����
 
    % ������������
    a0=-1i*2*pi*arrayEleDistance/lambda.*(0:arrayElementNumber-1);
    A0=exp(a0.'*sin(dec2rad(realDOA)));

    %��������źŲ���
    x=A0*signalsOutput+noiseMatrix;

    %�������R����
    R=(x*x')./sampleNumber;

    % ��R���������ֽ⣬�õ�U��Un��Us
    % ��������ʹ�ø���Ч��svd()������
    % ��ѧ�Ͽ���֤������ʹ��svd()��eig()��Ч
    [U,~]=svd(R);
    Un=U(:,signalNumber+1:arrayElementNumber);

    % ��ɨ�辫��ΪL
    L=500;
	% ��theta��-90�ȵ�90��֮��仯
    theta=-90:180/L:90-180/L;
    % ���ݹ�ʽ����P(\theta)
    a=exp(-1i.*2.*pi.*arrayEleDistance./lambda*(0:arrayElementNumber-1)'*sin(dec2rad(theta)));
    P=diag(1./(a'*(Un*Un')*a));

    % ��ͼ
    plot(theta,10*log10(abs(P/max(P))));grid;
    axis([-90,90,-20,0]);
    xlabel('\theta/degree');ylabel('P(\theta)/dB');
    title('P(\theta)');
end