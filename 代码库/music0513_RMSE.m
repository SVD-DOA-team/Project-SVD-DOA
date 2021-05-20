% ���ű�ʵ����ȷ���źŵ�ָ��������ι���ȡƽ����
K=5;
M=7;

%���ź�Ƶ��Ϊ13.56MHz
f=13.56*10^6;
c=3*10^8;
lmd=c/f;
d=lmd/2;
N=256;
t=0:N-1;
Ps=0;

% s=Ps*randn(K,N)+Ps*randn(K,N)*1i;
% ������N�������������±ض����þ����������ڵ������źš�
% ��������Ŀ���ж���һ�������źŵ������ʹ��ʵ������ظ������Ǳ��ⲻ��һ�����ڵ������źŷ�������й©�����⡣
% s=wgn(K,N,Ps);
s = zeros(K,N);
for index=1:K
    s(index,:)=10^(Ps/10)*sin(2*pi*index/(N-1)*t);
    % ͨ�����ַ�ʽ���ɵ��źţ���Ƶ�����ź������볤��(N-1)�ĺ���
    % ��������������֤���������ɵĺ����ܹ��ܺõ���N���ڱ�����������
    % �Ӷ������˹���й©
    % disp(['s',num2str(index),'=',num2str(s(index,N))])
end

% ��ʵ�Ĳ���� theta0
% Ϊ�����Զ�������������������ΪΪn*K*10��
theta0=pi/18.*(1:K);

maxWEqualsToIError = [];
maxWhenWLSIsUsedError = [];

% ���ڷ������ֹ��Ʒ���ƽ������ͳ������
maxavgErrorWequI = [];
maxavgErrorWLSIsUsed = [];
minavgErrorWequI = [];
minavgErrorWLSIsUsed = [];
RMSEWequI = [];
RMSEWLSIsUsed = [];


for SNR=0:55
    errorWequI100loops = [];
    errorWLS100loops = [];
    % ���ڼ������ֹ��Ʒ����Ĳв�ƽ����
    errorSqaureSumWequI = 0;
    errorSqaureSumWLSIsUsed = 0;
    loopTimes = 300;
    for loopCounter = 1:loopTimes
        % ��������ΪPn,��λdB
        Pn=Ps-SNR;

        % ��������
        n=wgn(M,N,Pn);

        %������������A0
        a0=1i*2*pi*d/lmd.*(0:M-1);
        A=exp(a0.'*sin(theta0));
        x=A*s+n;%��������źŲ���

        %�������R����
        R=(x*x')./N;
        Q=QMatrixGenerate(M);  %���� MxM ���ϱ任���� Q
        C=real(Q'*R*Q); % �ϱ任ȡʵ����MxM
        [U,sigma,V]=svd(C);% �� C �� SVD ׼��ȡ�ӿռ䣬MxM
        [~,LambdaR,~]=svd(R); % �� R ֱ���� SVD����������������ǲ����ģ�������
        LambdaR=LambdaR(1:K,1:K); % ȡǰ KxK �Ŀ飬�����������׵�
        Us=U(:,1:K);%MxK
        UC=Q*Us;%MxK
        D=zeros(K,(M-K)*K);
        f=zeros(1,(M-K)*K);
        for j=1:K
            fk=-UC(K+1:M,j);
            Dk=zeros(M-K,K);
            for i=1:K
                Dk(:,i)=UC(K-i+1:M-i,j);
            end
            D(:,(j-1)*(M-K)+1:j*(M-K))=Dk';
            f(:,(j-1)*(M-K)+1:j*(M-K))=fk';
        end
        D=D';%(M-K)*KxK
        f=f';%(M-K)*KxK

        % �˴���W=I������DOA
        W=eye((M-K)*K);
        c=(D'*W*D)\(D'*W*f);%Kx1
        c=[1;c];%(K+1)x1
        z=roots(c);%Kx1
        B=toeplitz([c(K+1),zeros(1,M-K-1)].',[(flip(c)).',zeros(1,M-K-1)]);
        W=kron(LambdaR.*LambdaR,inv(B*B'));
        u=log(z)/1i;%Kx1
        theta1=sort(asin(lmd/(2*pi*d)*u));%Kx1
        theta1=abs(real(theta1'));%1xK
        errorWhenWEqualsToI=abs(theta1-theta0)*180/pi;
        errorWequI100loops = [errorWequI100loops;errorWhenWEqualsToI];
        % ��¼�в��ƽ����
        errorSqaureSumWequI = errorSqaureSumWequI + errorWhenWEqualsToI.^2;
        
        % ������ʱ�����ĵ�����ÿһ���Ľ�����ʽ���һ�к͵�115��117�С�120��ע��
        % �����Ҫ�������Ľ��������һ�к͵�115�С�117�С�120�лָ�����
        errorWhenWLSIsUsed=[];
        Iter_num=6;
        % ��W��������ֵ
        W=eye((M-K)*K);
        for i=1:Iter_num
            c=(D'*W*D)\(D'*W*f);%Kx1
            c=[1;c];%(K+1)x1
            z=roots(c);%Kx1
            B=toeplitz([c(K+1),zeros(1,M-K-1)].',[(flip(c)).',zeros(1,M-K-1)]);
            W=kron(LambdaR.*LambdaR,inv(B*B'));
            u=log(z)/1i;%Kx1
            theta=sort(asin(lmd/(2*pi*d)*u));%Kx1
            theta=abs(real(theta'));%1xK
            errorWhenWLSIsUsed=[errorWhenWLSIsUsed;abs(theta-theta0)*180/pi];
        end
        errorWLS100loops = [errorWLS100loops;errorWhenWLSIsUsed(Iter_num,:)];
        % errorWLS100loops = [errorWLS100loops;abs(theta-theta0)*180/pi];
        % ��¼�в��ƽ����
        errorSqaureSumWLSIsUsed = errorSqaureSumWLSIsUsed + errorWhenWLSIsUsed(Iter_num,:).^2;
        % errorSqaureSumWLSIsUsed = errorSqaureSumWLSIsUsed + (abs(theta-theta0)*180/pi).^2;
    end
    index=1:K;
    averageErrorWequaltoI = mean(errorWequI100loops,1);
    averageErrorWLSIsUsed = mean(errorWLS100loops,1);
%     figure
%     plot(index,averageErrorWequaltoI,'b',index,averageErrorWLSIsUsed,'r')
%     title(['��SNRΪ',num2str(SNR),'ʱ���������ͼ']);
%     xlabel(['��������N=',num2str(loopTimes)])
%     ylabel('����С')
%     legend('WȡIʱ�����','WLS������W�����')
    maxavgErrorWequI = [maxavgErrorWequI max(averageErrorWequaltoI)];
    minavgErrorWequI = [minavgErrorWequI min(averageErrorWequaltoI)];
    
    maxavgErrorWLSIsUsed = [maxavgErrorWLSIsUsed max(averageErrorWLSIsUsed)];
    minavgErrorWLSIsUsed = [minavgErrorWLSIsUsed min(averageErrorWLSIsUsed)];
    
    % ��ÿһ��SNR��������������ʽ��¼(��1�Ĺ���RMSE, ��1�Ĺ���RMSE, ... ��5�Ĺ���RMSE)
    RMSEWequI = [RMSEWequI;sqrt(errorSqaureSumWequI/(loopTimes-1))];
    RMSEWLSIsUsed = [RMSEWLSIsUsed;sqrt(errorSqaureSumWLSIsUsed/(loopTimes-1))];
end
SNR=0:55;
% figure
% plot(SNR,maxavgErrorWequI,'m',SNR,maxavgErrorWLSIsUsed,'r',SNR,minavgErrorWequI,'b',SNR,minavgErrorWLSIsUsed,'k')
% title('��ͬSNR����������DOA���Ʒ�����ƽ�����')
% xlabel('SNR')
% ylabel('����С')
% legend('WȡIʱ��������','WLS������W��������','WȡIʱ����С���','WLS������W����С���')
% figure
% plot(SNR,maxavgErrorWequI,'m',SNR,maxavgErrorWLSIsUsed,'r')
% title('��ͬSNR����������DOA���Ʒ��������ƽ�����')
% xlabel('SNR')
% ylabel('����С')
% legend('WȡIʱ��������','WLS������W��������')
% figure
% plot(SNR,minavgErrorWequI,'b',SNR,minavgErrorWLSIsUsed,'k')
% title('��ͬSNR����������DOA���Ʒ�������Сƽ�����')
% xlabel('SNR')
% ylabel('����С')
% legend('WȡIʱ����С���','WLS������W����С���')
CRB = sqrt(loopTimes./10.^((SNR)./10));
for index = 1:K
    figure('Name',['����DOA���Ʒ������ƽ�',num2str(index),'�ľ��������'],'NumberTitle','off');
    plot(SNR,RMSEWequI(:,index),'m',SNR,RMSEWLSIsUsed(:,index),'b')
    title(['����DOA���Ʒ������ƽ�',num2str(index),'�ľ��������'])
    xlabel('SNR')
    ylabel('����С(RMSE)')
    legend('WȡIʱ��RMSE','WLS������W��RMSE')
end