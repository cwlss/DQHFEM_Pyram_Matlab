function [X,Y,Z,dX,dY,dZ,dJ]=MapCurTetra(r,s,t,x,y,z,number)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
if nargin<7
    number.edge=[0,0,0,0,0,0];
    number.face={0,0,0,0};
    number.H=0;
    x=x(1:4);
    y=y(1:4);
    z=z(1:4);
end

nr=length(r);
% Ne=sum(number.edge);
% Nf=0;
% for i=1:length(number.face)
%     Nf=Nf+number.face{i}(1)*number.face{i}(2);
% end
% N=Ne+Nf+8;
% if length(x)~=N;
%     error('x and number is not competitive')
% end

% x=MapRight{1};y=MapRight{2};z=MapRight{3};
if nargout==3
    G=BaseTetra(r,s,t,number);
    T=InterpTetra(number,2);
    GG=G*T;
    X=GG*x;
    Y=GG*y;
    Z=GG*z;
else
    [G,Gr,Gs,Gt]=BaseTetra(r,s,t,number );
    T=InterpTetra(number,2);
    [GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T);
    X=GG*x;
    Y=GG*y;
    Z=GG*z;
    dxdr=GR*x;
    dxds=GS*x;
    dxdt=GT*x;
    dX={dxdr,dxds,dxdt};
    dydr=GR*y;
    dyds=GS*y;
    dydt=GT*y;
    dY={dydr,dyds,dydt};
    dzdr=GR*z;
    dzds=GS*z;
    dzdt=GT*z;
    dZ={dzdr,dzds,dzdt};
    dJ=zeros(nr,1);
    for i=1:nr
        J=[dxdr(i),dydr(i),dzdr(i)
            dxds(i),dyds(i),dzds(i)
            dxdt(i),dydt(i),dzdt(i)];
        dJ(i)=det(J);
    end
end