function [Gx,Gy,Gz]=MapTrans(Gr,Gs,Gt,dX,dY,dZ,dJ)
% 链式法则导致的坐标变换，参见专著7.6节
num=length(dJ);
dxdr=dX{1};dxds=dX{2};dxdt=dX{3};
dydr=dY{1};dyds=dY{2};dydt=dY{3};
dzdr=dZ{1};dzds=dZ{2};dzdt=dZ{3};

Gx=zeros(size(Gs));Gy=Gx;Gz=Gx;
for i=1:num
jx=1/dJ(i)*[dyds(i)*dzdt(i)-dydt(i)*dzds(i),dydt(i)*dzdr(i)-dydr(i)*dzdt(i),dydr(i)*dzds(i)-dyds(i)*dzdr(i)];
jy=1/dJ(i)*[dxdt(i)*dzds(i)-dxds(i)*dzdt(i),dxdr(i)*dzdt(i)-dxdt(i)*dzdr(i),dxds(i)*dzdr(i)-dxdr(i)*dzds(i)];
jz=1/dJ(i)*[dxds(i)*dydt(i)-dxdt(i)*dyds(i),dxdt(i)*dydr(i)-dxdr(i)*dydt(i),dxdr(i)*dyds(i)-dxds(i)*dydr(i)];
Gx(i,:)=jx*[Gr(i,:);Gs(i,:);Gt(i,:)];
Gy(i,:)=jy*[Gr(i,:);Gs(i,:);Gt(i,:)];
Gz(i,:)=jz*[Gr(i,:);Gs(i,:);Gt(i,:)];
end
end
