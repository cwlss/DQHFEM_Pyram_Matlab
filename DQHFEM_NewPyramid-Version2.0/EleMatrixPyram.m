function  [K,M] =EleMatrixPyram(ElementData,E,mum,rhom)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明

number=ElementData.Number;
eleX=ElementData.MapRight;
if length(eleX{1})==5
MapNumber.edge=[0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],0,0,0,0};
MapNumber.H=0;
else
    MapNumber=ElementData.MapNumber;
end
%% 给点积分点
nq=IntegOrderPyram(number);
IP=IntegralPntPyram(nq);
%% 求取积分点处的形函数值
[G,Gr,Gs,Gt]=BasePyram(IP.Ri,IP.Si,IP.Ti,number);
%% 将边界上的形函数转化为插值形式
T=InterpPyram(number);
[GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T);
%% 将形函数关于自然坐标r,s,t的导数值映射为关于真实坐标x,y,z的导数值
x=eleX{1};y=eleX{2};z=eleX{3};
[~,~,~,dX,dY,dZ,dJ]=MapCurPyram(IP.Ri,IP.Si,IP.Ti,x,y,z,MapNumber);
[Gx,Gy,Gz]=MapTrans(GR,GS,GT,dX,dY,dZ,dJ);
%% 求刚度矩阵和质量矩阵
D=ElasticMatrix(E,mum);%弹性矩阵
B=StrainMatrix(Gx,Gy,Gz);%应变矩阵
K=StiffnessMatrix(B,D,dJ,IP.Ci);%刚度矩阵
if nargout>1
M=MassMatrix(GG,dJ,IP.Ci,rhom);%质量矩阵
end
end

function nq=IntegOrderPyram(number)
pz1=max(number.edge(1:4))+2;
pz2=max(number.edge(5:8))+2;
pz3=sum(number.face{1})+2;
pz4=max([number.face{2:5}])+3;
pz5=number.H+4;
pz=max([pz1,pz2,pz3,pz4,pz5]);

px1=max(number.edge([1,3]))+2;
px2=max(number.face{1}(1))+3;
px3=max([number.face{[2,4]}])+3;
px4=number.H+3;
px=max([px1,px2,px3,px4]);

py1=max(number.edge([2,4]))+2;
py2=max(number.face{1}(2))+3;
py3=max([number.face{[1,3]}])+3;
py4=number.H+3;
py=max([py1,py2,py3,py4]);
nq=[px+1,py+1,pz+1];
end
    