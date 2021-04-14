
%music����ʵ��
%���貨�﷽���\theta0;
%����ֻ�м��԰�����N(t)
%���ź�Ƶ��Ϊ13.56MHz
%��Ϊ����Ԫ�������ʱ��������Ƚ��ͣ�����ѡ����Ԫ���dΪ\lamda/2
%�ź�Դ��K����Ԫ��M��
K=2;
M=4;
f=13.56*10^6;
c=3*10^8;
lmd=c/f;
d=lmd/2;
N=120;
t=0:N-1;
%SNR��λdB
SNR=1.5;
A=10;
Ps=A^2/2;
w=[1/3*pi;1/6*pi];
%theta0Ϊn*60/K��
theta0=60./(1:K);
%��Դ�������ź�
s=10*cos(w*t);
%��������ΪN0
N0=Ps/10^SNR;
n=N0*randn(M,N)+N0*randn(M,N)*1i;

%������������A0
a0=-1i*2*pi*d/lmd.*(0:M-1);
A0=exp(a0.'*sin(theta0/180*pi));

x=A0*s+n;%��������źŲ���

%�������R����
R=(x*x')./N;

%��R����svd�ֽ��U��Un��Us
[U,sigma]=svd(R);
Us=U(:,1:K);
Un=U(:,K+1:M);

%��theta��-90��-90��֮��仯
%��ɨ�辫��ΪL
%���ݹ�ʽ����P(\theta)
L=10000;
theta=-90:180/L:90-180/L;
a=exp(-1i.*2.*pi.*d./lmd*(0:M-1)'*sin(theta/180*pi));
P=diag(1./(a'*(Un*Un')*a));



%���ȷ���
%����
F = 10*log10(abs(P/max(P)));
[pks,locs]=findpeaks(F,theta); %pks��ֵ
% num=length(locs);%�жϵó��ķ��ֵλ�õ�����
% j = 1;
%�˳�����
THETA=find(abs(pks)<=30);%ȷ����Ӧ��psk������
PSK = pks(find(abs(pks)<=30))
num=length(THETA);
theta1 = zeros(1,num);%����DOA���ƵĽǶȴ洢����
% PSK =zeros(1,num);%�˳����psk
%����洢
for i = 1: num
    theta1(i)=locs(THETA(i));
%     PSK(i)=F(theta1(i));
end
%���ȹ���

pks
theta1

%��ͼ
% % figure
% % plot(theta,F,theta0,PSK,'x');grid;
% % axis([-90,90,-20,0]);
% % xlabel("\theta/degree");ylabel("P(\theta)/dB");
% % title("P(\theta)");
%����M����߷�ֵ��Ԥ�⾫�ȣ���������Ȼ��С�����ȣ����Ԥ�⾫�ȣ�����δ���ж�����ͼ������
%������A0��s��a�����ɷ�ʽ������˳�������Ч��