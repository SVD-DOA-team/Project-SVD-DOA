function Q=QMatrixGenerate(M)
if(mod(M,2)==0)
    l=M/2;
    JL=flip(eye(l));
    IL=eye(l);
    Q=sqrt(2)/2*[[IL,1i*IL];[JL,-1i*JL]];
else
    l=(M-1)/2;
    JL=flip(eye(l));
    IL=eye(l);
    Z=zeros(l,1);
    Q=sqrt(2)/2*[[IL,Z,1i*IL];[Z',sqrt(2),Z'];[JL,Z,-1i*JL]];
end
end