function [X,Y,Z,dX,dY,dZ,dJ]=MapCurHex(r,s,t,x,y,z,number)
% 缺省调用1：[X,Y,Z,dX,dY,dZ,dJ]=MapCurHex(r,s,t,x,y,z) 为直边映射
% 缺省调用2：[X,Y,Z]=MapCurHex(r,s,t,x,y,z，number) 
% 缺省调用3：[X,Y,Z]=MapCurHex(r,s,t,x,y,z) 直边映射

% r，s，t：参数点坐标（积分点）
% x,y,z映射权矢量（线性映射为单元的顶点几何坐标值，曲边映射为边界映射节点的几何坐标值）
%number：结构体：储存映射的阶次（与各边界映射节点数量相一致）
% 注意[这里的number和基函数中的number，虽然字段相同，但表示的意义不同，不再是单元位移逼近的阶次，而是单元形状几何逼近的阶次]

% X,Y,Z：几何域坐标（积分点在几何域中的坐标）
% dX,dY,dZ参考专著7.6节
% dJ 雅各比行列式的值
if nargin<7
    number.edge=[0,0,0,0,0,0,0,0,0,0,0,0];
    number.face={[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]};
    number.H=[0,0,0];
    x=x(1:8);
    y=y(1:8);
    z=z(1:8);
end
nr=length(r);

if nargout==3
    G=BaseHex(r,s,t,number);
    T=InterpHex(number,2);
    GG=G*T;
    X=GG*x;
    Y=GG*y;
    Z=GG*z;
else
    [G,Gr,Gs,Gt]=BaseHex(r,s,t,number );
    T=InterpHex(number,2);
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