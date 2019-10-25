function  [K,M] =EleMatrixHex(ElementData,E,mum,rhom)
%UNTITLED2 �˴���ʾ�йش˺�����ժҪ
% ��Ԫ����K M

%   �˴���ʾ��ϸ˵��
number=ElementData.Number;
eleX=ElementData.MapRight;
if length(eleX{1})==8
MapNumber.edge=[0,0,0,0,0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],[0,0],[0,0],[0,0],[0,0],[0,0]};
MapNumber.H=[0,0,0];
else
    MapNumber=ElementData.MapNumber;
end
    
%% ������ֵ�
nq=OrderIntegral(number);
IP=IntegralPntHex(nq);
%% ��ȡ���ֵ㴦���κ���ֵ
[G,Gr,Gs,Gt]=BaseHex(IP.Ri,IP.Si,IP.Ti,number);
%% ���߽��ϵ��κ���ת��Ϊ��ֵ��ʽ
T=InterpHex(number);
[GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T);
%% ���κ���������Ȼ����r,s,t�ĵ���ֵӳ��Ϊ������ʵ����x,y,z�ĵ���ֵ
x=eleX{1};y=eleX{2};z=eleX{3};
[~,~,~,dX,dY,dZ,dJ]=MapCurHex(IP.Ri,IP.Si,IP.Ti,x,y,z,MapNumber);
[Gx,Gy,Gz]=MapTrans(GR,GS,GT,dX,dY,dZ,dJ);
%% ��նȾ������������
D=ElasticMatrix(E,mum);%���Ծ���
B=StrainMatrix(Gx,Gy,Gz);%Ӧ�����
K=StiffnessMatrix(B,D,dJ,IP.Ci);%�նȾ���
if nargout>1
M=MassMatrix(GG,dJ,IP.Ci,rhom);%��������
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
    