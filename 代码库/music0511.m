% 本脚本实现了确定信号的指定次数多次估计取平均。
K=5;
M=7;

%设信号频率为13.56MHz
f=13.56*10^6;
c=3*10^8;
lmd=c/f;
d=lmd/2;
N=256;
t = 0:N-1;
Ps=0;

% s=Ps*randn(K,N)+Ps*randn(K,N)*1i;
% 改用在N点采样点数情况下必定正好经过完整周期的正弦信号。
% 这样做的目的有二，一是消除信号的随机性使得实验可以重复，二是避免不满一个周期的正弦信号发生功率泄漏的问题。
% s=wgn(K,N,Ps);
s = zeros(K,N);
for index=1:K
    s(index,:)=10^(Ps/10)*sin(2*pi*index/(N-1)*t);
    % 通过这种方式生成的信号，其频率是信号索引与长度(N-1)的函数
    % 下面这条语句可以证明这样生成的函数能够很好地在N点内保持完整周期
    % 从而避免了功率泄漏
    % disp(['s',num2str(index),'=',num2str(s(index,N))])
end

% 真实的波达角 theta0
% 为便于自动生成起见，这里简单设置为为n*K*10度
theta0=pi/18.*(1:K);

maxWEqualsToIError = [];
maxWhenWLSIsUsedError = [];
for SNR=[2:2:8 10:10:60]
    errorWequI100loops = [];
    errorWLS100loops = [];
    for loopTimes = 1:1000
        % 噪声功率为Pn,单位dB
        Pn=Ps-SNR;

        % 生成噪声
        n=wgn(M,N,Pn);

        %生成阵列流形A0
        a0=1i*2*pi*d/lmd.*(0:M-1);
        A=exp(a0.'*sin(theta0));
        x=A*s+n;%阵列输出信号采样

        %求出采样R矩阵
        R=(x*x')./N;
        Q=QMatrixGenerate(M);  %生成 MxM 的酉变换矩阵 Q
        C=real(Q'*R*Q); % 酉变换取实部，MxM
        [U,sigma,V]=svd(C);% 对 C 作 SVD 准备取子空间，MxM
        [~,LambdaR,~]=svd(R); % 对 R 直接作 SVD，左右奇异矩阵我们不关心，故抛弃
        LambdaR=LambdaR(1:K,1:K); % 取前 KxK 的块，将噪声部分抛掉
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

        % 此处用W=I来估计DOA
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

        errorWhenWLSIsUsed=[];
        Iter_num=50;
        % 给W赋迭代初值
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
    end
    index=1:K;
    figure
    averageErrorWequaltoI = mean(errorWequI100loops,1);
    averageErrorWLSIsUsed = mean(errorWLS100loops,1);
    plot(index,averageErrorWequaltoI,'b',index,averageErrorWLSIsUsed,'r')
    title(['当SNR为',num2str(SNR),'时的误差折线图']);
    xlabel(['迭代次数N=',num2str(loopTimes)])
    ylabel('误差大小')
    legend('W取I时的误差','WLS迭代求W的误差')
end
