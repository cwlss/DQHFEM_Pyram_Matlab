function  [K,M] =EleMatrixHex(ElementData,E,mum,rhom)
%UNTITLED2 此处显示有关此函数的摘要
% 单元矩阵K M

%   此处显示详细说明
number=ElementData.Number;
eleX=ElementData.MapRight;
if length(eleX{1})==8
MapNumber.edge=[0,0,0,0,0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]};
MapNumber.H=[0,0,0];
else
    MapNumber=ElementData.MapNumber;
end
    
%% 给点积分点
nq=OrderIntegral(number);
IP=IntegralPntHex(nq);
%% 求取积分点处的形函数值
[G,Gr,Gs,Gt]=BaseHex(IP.Ri,IP.Si,IP.Ti,number);
%% 将边界上的形函数转化为插值形式
T=InterpHex(number);
[GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T);
%% 将形函数关于自然坐标r,s,t的导数值映射为关于真实坐标x,y,z的导数值
x=eleX{1};y=eleX{2};z=eleX{3};
[~,~,~,dX,dY,dZ,dJ]=MapCurHex(IP.Ri,IP.Si,IP.Ti,x,y,z,MapNumber);
[Gx,Gy,Gz]=MapTrans(GR,GS,GT,dX,dY,dZ,dJ);
%% 求刚度矩阵和质量矩阵
D=ElasticMatrix(E,mum);%弹性矩阵
B=StrainMatrix(Gx,Gy,Gz);%应变矩阵
K=StiffnessMatrix(B,D,dJ,IP.Ci);%刚度矩阵
if nargout>1
M=MassMatrix(GG,dJ,IP.Ci,rhom);%质量矩阵
end
end

function nq=OrderIntegral(number)
MAX=max(number.edge);
for i=1:length(number.face)
    Max=max(number.face{i});
    if MAX<Max
        MAX=Max;
    end
end
if MAX<max(number.H)
    MAX=max(number.H);
end
dq=4;
nq=[MAX+dq,MAX+dq,MAX+dq];
end
    