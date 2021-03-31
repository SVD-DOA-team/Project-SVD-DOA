clc;clear; close all;
%% ���ɲ��Ի���
%music����ʵ��
%���貨�﷽���\theta0;
%����ֻ�м��԰�����N(t)
%���ź�Ƶ��Ϊ13.56MHz
%��Ϊ����Ԫ�������ʱ��������Ƚ��ͣ�����ѡ����Ԫ���dΪ\lamda
%�ź�Դ��K����Ԫ��M��
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
    s0(index,:)=10*cos(w(index).*t);%��Դ�������ź�
end
%������ֵ���ΪN0
% ���������ķ���
% N0=10;
N0=10;
% ȡһ���̶�����
n=N0*randn(M,N)+N0*randn(M,N)*1i;
A0=zeros(M,K);
for index=1:M
    for h=1:K
        A0(index,h)=exp((-1i*2*pi*d/lmd*sin(theta0(h)))*(index-1));
        %disp(['A0(',num2str(index),',',num2str(h),')=',num2str(A0(index,h))])
    end
end
%% ������������źŽ��
% % �� s0 ���˳�һ������Ȥ�� s
% for i=1:K
%    s=[zeros(i-1,N);s0(i,:);zeros(K-i,N)];
%    % ��ʵ�ֵĴ���
% x=A0*s+n;%��������źŲ���
% 
% %�������R����
% R=(x*x')./N;
% 
% %��R����svd�ֽ��U��Un��Us
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
% %��theta��0��pi֮��仯
% %��ɨ�辫��ΪL
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
%% �����֣��������
x=A0*s0+n;%��������źŲ���

%�������R����
R=(x*x')./N;

%��R����svd�ֽ��U��Un��Us
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

%��theta��0��pi֮��仯
%��ɨ�辫��ΪL
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
