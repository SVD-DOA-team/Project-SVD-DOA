%music基础实现
%干扰只有加性白噪声N(t)
%因为当阵元间距扩大时，波束宽度降低，所以选择阵元间距d为\lambda/2
%信号源数K；阵元数M；
K=2;
M=5;

%设信号频率为13.56MHz
f=13.56*10^6;

% 光速c，常数
c=3*10^8;

% 计算得到波长
lmd=c/f;

% 阵元间距
% 这个值的选择会影响到DOA估计，可能与空域采样的奈奎斯特定理有关。
d=lmd/2;

% 时域采样点数/快拍数
N=120;

% SNR，单位dB
SNR=5;

% 信号功率,单位dB
Ps=0;

% s=Ps*randn(K,N)+Ps*randn(K,N)*1i;
% 因为DOA估计效果仅与信源功率有关，同时减少正弦信号采样不为周期整数倍时的功率误差，故用零均值高斯白噪声来做复信源
s=wgn(K,N,Ps);
% 用白噪声生成公式优化

% 真实的波达角 theta0
% 为便于自动生成起见，这里简单设置为为n*K*10度
theta0=10.*(1:K);

% 噪声功率为Pn,单位dB
Pn=Ps-SNR;

% 生成噪声
% n=N0*randn(M,N)+N0*randn(M,N)*1i;
n=wgn(M,N,Pn);

%生成阵列流形A0
a0=-1i*2*pi*d/lmd.*(0:M-1);
A=exp(a0.'*sin(theta0/180*pi));

x=A*s+n;%阵列输出信号采样

%求出采样R矩阵
% R=(x*x')./N;

Rs=10^(Ps/10)*eye(K,K);
R=A*Rs*A'+10^(Pn/10)*eye(M,M);
%完美R矩阵

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

%画图
plot(theta,10*log10(abs(P/max(P))));grid;
axis([-90,90,-60,0]);
xlabel("\theta/degree");ylabel("P(\theta)/dB");
title("P(\theta)");
