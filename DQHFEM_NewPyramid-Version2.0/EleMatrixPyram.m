function  [K,M] =EleMatrixPyram(ElementData,E,mum,rhom)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��

number=ElementData.Number;
eleX=ElementData.MapRight;
if length(eleX{1})==5
MapNumber.edge=[0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],0,0,0,0};
MapNumber.H=0;
else
    MapNumber=ElementData.MapNumber;
end
%% ������ֵ�
nq=IntegOrderPyram(number);
IP=IntegralPntPyram(nq);
%% ��ȡ���ֵ㴦���κ���ֵ
[G,Gr,Gs,Gt]=BasePyram(IP.Ri,IP.Si,IP.Ti,number);
%% ���߽��ϵ��κ���ת��Ϊ��ֵ��ʽ
T=InterpPyram(number);
[GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T);
%% ���κ���������Ȼ����r,s,t�ĵ���ֵӳ��Ϊ������ʵ����x,y,z�ĵ���ֵ
x=eleX{1};y=eleX{2};z=eleX{3};
[~,~,~,dX,dY,dZ,dJ]=MapCurPyram(IP.Ri,IP.Si,IP.Ti,x,y,z,MapNumber);
[Gx,Gy,Gz]=MapTrans(GR,GS,GT,dX,dY,dZ,dJ);
%% ��նȾ������������
D=ElasticMatrix(E,mum);%���Ծ���
B=StrainMatrix(Gx,Gy,Gz);%Ӧ�����
K=StiffnessMatrix(B,D,dJ,IP.Ci);%�նȾ���
if nargout>1
M=MassMatrix(GG,dJ,IP.Ci,rhom);%��������
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
    