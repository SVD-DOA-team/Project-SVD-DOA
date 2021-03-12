% clc;clear;close all;
% errU=[];
% errS=[];
% errV=[];
% for i=1:1
%    a=fix(20*rand)+2;
%    b=fix(20*rand)+2;
%    mat=rand(a,b);
%    [stdU,stdS,stdV]=svd(mat);
%    [myU,myS,myV]=mysvd(mat);
% %    errU=[errU,max(abs(myU-stdU))];
% %    errS=[errS,max(abs(myS-stdS))];
% %    errV=[errV,max(abs(myV-stdV))];
%     
%    errU=max(max(abs(myU-stdU)));
%    errS=max(max(abs(myS-stdS)));
%    errV=max(max(abs(myV-stdV)));
%    if max(errU,errS,errV)>1e-6
%        pause(inf);
%    end
% end
% % max(errU)
% % max(errS)
% % max(errV)





clc;clear;close all;
errU=[];
errS=[];
errV=[];
for i=1:1000
   a=fix(20*rand)+2;
   b=fix(20*rand)+2;
   mat=rand(a,b);
   [stdU,stdS,stdV]=svd(mat);
   [myU,myS,myV]=mysvd(mat);
   errU=[errU,max(abs(myU-stdU))];
   errS=[errS,max(abs(myS-stdS))];
   errV=[errV,max(abs(myV-stdV))];
end
errU
errS
errV