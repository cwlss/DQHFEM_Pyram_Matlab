function [X,Y,Z,dX,dY,dZ,dJ]=MapCurHex(r,s,t,x,y,z,number)
% ȱʡ����1��[X,Y,Z,dX,dY,dZ,dJ]=MapCurHex(r,s,t,x,y,z) Ϊֱ��ӳ��
% ȱʡ����2��[X,Y,Z]=MapCurHex(r,s,t,x,y,z��number) 
% ȱʡ����3��[X,Y,Z]=MapCurHex(r,s,t,x,y,z) ֱ��ӳ��

% r��s��t�����������꣨���ֵ㣩
% x,y,zӳ��Ȩʸ��������ӳ��Ϊ��Ԫ�Ķ��㼸������ֵ������ӳ��Ϊ�߽�ӳ��ڵ�ļ�������ֵ��
%number���ṹ�壺����ӳ��Ľ״Σ�����߽�ӳ��ڵ�������һ�£�
% ע��[�����number�ͻ������е�number����Ȼ�ֶ���ͬ������ʾ�����岻ͬ�������ǵ�Ԫλ�Ʊƽ��Ľ״Σ����ǵ�Ԫ��״���αƽ��Ľ״�]

% X,Y,Z�����������꣨���ֵ��ڼ������е����꣩
% dX,dY,dZ�ο�ר��7.6��
% dJ �Ÿ�������ʽ��ֵ
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