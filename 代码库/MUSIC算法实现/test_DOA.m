clc;clear; close all;
%% 生成测试环境
%music基础实现
%假设波达方向角\theta0;
%干扰只有加性白噪声N(t)
%设信号频率为13.56MHz
%因为当阵元间距扩大时，波束宽度降低，所以选择阵元间距d为\lamda
%信号源数K；阵元数M；
K=6;
M=18;
f=13.56*10^6;
c=3*10^8;
lmd=c/f;
d=lmd/2;
N=1200;
t=0:N-1;
w=(pi/3)./(1:K);
theta0=pi/12.*(1:K);
s0=zeros(K,N);
for index=1:K
    s0(index,:)=10*cos(w(index).*t);%信源发出的信号
end
%噪声幅值最大为N0
% 降低噪声的幅度
% N0=10;
N0=10;
% 取一个固定噪声
n=N0*randn(M,N)+N0*randn(M,N)*1i;
A0=zeros(M,K);
for index=1:M
    for h=1:K
        A0(index,h)=exp((-1i*2*pi*d/lmd*sin(theta0(h)))*(index-1));
        %disp(['A0(',num2str(index),',',num2str(h),')=',num2str(A0(index,h))])
    end
end
%% 独立输出六个信号结果
% % 从 s0 中滤出一个感兴趣的 s
% for i=1:K
%    s=[zeros(i-1,N);s0(i,:);zeros(K-i,N)];
%    % 待实现的代码
% x=A0*s+n;%阵列输出信号采样
% 
% %求出采样R矩阵
% R=(x*x')./N;
% 
% %对R进行svd分解出U，Un，Us
% [U,sigma]=svd(R);
% Us=zeros(M,1);
% K0=10;
% for i=1:K0
%     if Us==0
%         Us=U(:,i);
%     else
%         Us=[Us,U(:,i)];
%     end
% end
% Un=zeros(M,1);
% for i=K0+1:M
%     if Un==0
%         Un=U(:,i);
%     else
%         Un=[Un,U(:,i)];
%     end
% end
% 
% %让theta在0：pi之间变化
% %设扫描精度为L
% L=10000;
% theta=0:pi/L:pi-pi/L;
% a=zeros(M,L);
% for index=1:M
%     a(index,:)=exp((-1i.*2.*pi.*d./lmd.*sin(theta))*(index-1));
% end
% P=diag(1./(a'*Un*Un'*a));
% figure;
% plot(theta/pi,10*log10(abs(P/max(P))));grid;
% axis([0,1,-30,-3])
% xlabel('\theta/\pi');ylabel('P(\theta)/dB');
% title('P(\theta)');
% end
%% 不区分，集体输出
x=A0*s0+n;%阵列输出信号采样

%求出采样R矩阵
R=(x*x')./N;

%对R进行svd分解出U，Un，Us
[U,sigma]=svd(R);
Us=zeros(M,1);
K0=10;
for i=1:K0
    if Us==0
        Us=U(:,i);
    else
        Us=[Us,U(:,i)];
    end
end
Un=zeros(M,1);
for i=K0+1:M
    if Un==0
        Un=U(:,i);
    else
        Un=[Un,U(:,i)];
    end
end

%让theta在0：pi之间变化
%设扫描精度为L
L=10000;
theta=0:pi/L:pi-pi/L;
a=zeros(M,L);
for index=1:M
    a(index,:)=exp((-1i.*2.*pi.*d./lmd.*sin(theta))*(index-1));
end
P=diag(1./(a'*Un*Un'*a));
plot(theta/pi,10*log10(abs(P/max(P))));grid;
axis([0,1,-30,-3])
xlabel('\theta/\pi');ylabel('P(\theta)/dB');
title('P(\theta)');
