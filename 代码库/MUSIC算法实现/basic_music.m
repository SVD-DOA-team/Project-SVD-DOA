%music基础实现
%假设波达方向角\theta0为60°
%干扰只有加性白噪声N(t)
%设信号频率为13.56MHz
%因为当阵元间距扩大时，波束宽度降低，所以选择阵元间距d为lambda
%信号源数K=1；阵元数M=2；
K=1;
M=2;
theta0=1/3*pi;
f=13.56*10^6;
c=3*10^8;
lambda=c/f;
d=lambda;
N=120;
t=0:N-1;
s=10*cos(2/15*pi.*t);%信源发出的信号
n=rand(2,N);
a0=zeros(2,1);
for i=1:2
    a0(i)=exp((-1i*2*pi*d/lambda*sin(theta0))*(i-1));
end
x=a0*s+n;%阵列输出信号采样

%求出采样R矩阵
R=(x*x')./N;

%对R进行svd分解出U，Un，Us
[U,sigma]=svd(R);
Us=zeros(M,1);
for i=1:K
    if Us==0
        Us=U(:,i);
    else
        Us=[Us,U(:,i)];
    end
end
Un=zeros(M,1);
for i=K+1:M
    if Un==0
        Un=U(:,i);
    else
        Un=[Un,U(:,i)];
    end
end

%让theta在0：pi之间变化
theta=0:0.0001*pi:pi-0.0001*pi;
a=zeros(2,10000);
for i=1:2
    a(i,:)=exp((-1i.*2.*pi.*d./lambda.*sin(theta))*(i-1));
end
P=diag(1./(a'*Un*Un'*a));
plot(theta/pi,10*log10(P/max(P)));grid;
axis([0,1,-10,0]);
xlabel('\theta/\pi');
ylabel('P(\theta)/dB');
title('P(\theta)');
