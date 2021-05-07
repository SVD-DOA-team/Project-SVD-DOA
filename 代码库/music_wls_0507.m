% MUSIC - 加权最小二乘第一版实现
% 目前的问题是程序会自发地在迭代2-4次后收敛，收敛位置的估计值不满足精度要求

%因为当阵元间距扩大时，波束宽度降低，所以选择阵元间距d为\lamda/2
%信号源数K；阵元数M；
K=5;
M=7;

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
SNR=60;
% 信号功率,单位dB
Ps=0;

% s=Ps*randn(K,N)+Ps*randn(K,N)*1i;
% 因为DOA估计效果仅与信源功率有关，同时减少正弦信号采样不为周期整数倍时的功率误差，故用零均值高斯白噪声来做复信源
s=wgn(K,N,Ps);
% 用白噪声生成公式优化

% 真实的波达角 theta0
% 为便于自动生成起见，这里简单设置为为n*K*10度
theta0=pi/18.*(1:K);

% 噪声功率为Pn,单位dB
Pn=Ps-SNR;

% 生成噪声
% n=N0*randn(M,N)+N0*randn(M,N)*1i;
n=wgn(M,N,Pn);
% n=zeros(M,N); %无噪情况下

%生成阵列流形A0
a0=1i*2*pi*d/lmd.*(0:M-1);
A=exp(a0.'*sin(theta0));

x=A*s+n;%阵列输出信号采样

%求出采样R矩阵
R=(x*x')./N;

% Rs=10^(Ps/10)*eye(K,K);
% R=A*Rs*A'+10^(Pn/10)*eye(M,M);
% %完美R矩阵
% 
% %对R进行svd分解出U，Un，Us
% [U,sigma]=svd(R);
% Us=U(:,1:K);
% Un=U(:,K+1:M);
% 
% %让theta在-90度-90度之间变化
% %设扫描精度为L
% %依据公式计算P(\theta)
% L=10000;
% theta=-90:180/L:90-180/L;
% a=exp(-1i.*2.*pi.*d./lmd*(0:M-1)'*sin(theta/180*pi));
% P=diag(1./(a'*(Un*Un')*a));
% 
% %画图
% plot(theta,10*log10(abs(P/max(P))));grid;
% axis([-90,90,-60,0]);
% xlabel("\theta/degree");ylabel("P(\theta)/dB");
% title("P(\theta)");

Q=QMatrixGenerate(M);  %MxM
C=real(Q'*R*Q); %MxM
[U,sigma,V]=svd(C);%MxM
[~,LambdaR,~]=svd(R);
LambdaR=LambdaR(1:K,1:K);
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

% 给W赋迭代初值
W=eye((M-K)*K);

% 初始化error
error=[];
for i=1:20
    c=(D'*W*D)\(D'*W*f);%Kx1
    c=[1;c];%(K+1)x1
    z=roots(c);%Kx1
    B=toeplitz([c(K+1),zeros(1,M-K-1)].',[(flip(c)).',zeros(1,M-K-1)]);
    W=kron(LambdaR.*LambdaR,inv(B*B'));
    u=log(z)/1i;%Kx1
	theta=sort(asin(lmd/(2*pi*d)*u));%Kx1
    theta=theta';%1xK
    error=[error;abs(theta-theta0)];
end
