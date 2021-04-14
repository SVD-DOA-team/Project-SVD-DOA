
%music基础实现
%假设波达方向角\theta0;
%干扰只有加性白噪声N(t)
%设信号频率为13.56MHz
%因为当阵元间距扩大时，波束宽度降低，所以选择阵元间距d为\lamda/2
%信号源数K；阵元数M；
K=2;
M=4;
f=13.56*10^6;
c=3*10^8;
lmd=c/f;
d=lmd/2;
N=120;
t=0:N-1;
%SNR单位dB
SNR=1.5;
A=10;
Ps=A^2/2;
w=[1/3*pi;1/6*pi];
%theta0为n*60/K度
theta0=60./(1:K);
%信源发出的信号
s=10*cos(w*t);
%噪声功率为N0
N0=Ps/10^SNR;
n=N0*randn(M,N)+N0*randn(M,N)*1i;

%生成阵列流形A0
a0=-1i*2*pi*d/lmd.*(0:M-1);
A0=exp(a0.'*sin(theta0/180*pi));

x=A0*s+n;%阵列输出信号采样

%求出采样R矩阵
R=(x*x')./N;

%对R进行svd分解出U，Un，Us
[U,sigma]=svd(R);
Us=U(:,1:K);
Un=U(:,K+1:M);

%让theta在-90度-90度之间变化
%设扫描精度为L
%依据公式计算P(\theta)
L=10000;
theta=-90:180/L:90-180/L;
a=exp(-1i.*2.*pi.*d./lmd*(0:M-1)'*sin(theta/180*pi));
P=diag(1./(a'*(Un*Un')*a));



%精度分析
%求峰峰
F = 10*log10(abs(P/max(P)));
[pks,locs]=findpeaks(F,theta); %pks峰值
% num=length(locs);%判断得出的峰峰值位置的数量
% j = 1;
%滤除干扰
THETA=find(abs(pks)<=30);%确定对应的psk的序列
PSK = pks(find(abs(pks)<=30))
num=length(THETA);
theta1 = zeros(1,num);%创建DOA估计的角度存储数组
% PSK =zeros(1,num);%滤除后的psk
%数组存储
for i = 1: num
    theta1(i)=locs(THETA(i));
%     PSK(i)=F(theta1(i));
end
%精度估计

pks
theta1

%画图
% % figure
% % plot(theta,F,theta0,PSK,'x');grid;
% % axis([-90,90,-20,0]);
% % xlabel("\theta/degree");ylabel("P(\theta)/dB");
% % title("P(\theta)");
%增大M会提高峰值与预测精度，增大信噪比会减小主瓣宽度，提高预测精度；但都未进行定量画图分析；
%更改了A0，s，a的生成方式，提高了程序运行效率