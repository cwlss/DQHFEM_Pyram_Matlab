%%2019.07.23 by Guomao
clear;clc
%% �������ͣ�ģ̬
% ����ֵ����㷨
solveway=2;%������ⷨ�������������ɱ߽�����
%% ǰ����ģ��
% ���ϲ���
E=195E9;
mum=0.3;
rhom=7722.3;
% rhom=100;
% �����ߴ�
a=100;b=100;h=100;%����
R=100;H=20;m=15;%Բ�壺mΪ����ӳ��ڵ�ĵ����������m����ȡ20���£�̫��ᵼ�²�ֵת������ӽ����죨���ھ���Ԫ�ؾ�С��1������ʽά��Խ��ֵ�ӽ���0��
% ��Ԫ�״�
for n=2
    %##################################################################
    % ElementDataԪ�����󣺳���Ϊ��Ԫ������Ԫ��Ϊ������Ԫ��Ӧ�ṹ�壬ÿ���ṹ�干�е��ֶ�����
    % 'Type'�ֶΣ���Ԫ������Ϣ
    % 'Number'�ֶΣ���Ԫ�߽�Ļ�������ֵ�ڵ���������״�n���
    % 'MapRight'�ֶ�:ӳ��Ȩϵ����Ϊ��Ԫӳ��ڵ�ļ������꣨ע����˳��Ҫ��
    % 'MapNumber'�ֶΣ���Ԫ�߽�ļ���ӳ��ڵ���������߽�ѡȡ��ӳ��ڵ���Ŀ��λ�����
    
% ���¾�Ϊģ��¼�����ǰ������̣���ȡ��ע�ͼ��ɸ����������м��㣬
% �����Ҫ�½���ģ�ͣ����½���һ��¼����򼴿ɣ�
% ����¼��������в���Ҫ�ģ��ɴӳ������ɾ��������Ӱ���������̵ļ���

    ElementData=huojian(n);%���
% ElementData=CirclePlate1(R,H,m,n);%Բ��ӳ�䲻��ȫ����
%     ElementData=CirPlate1(R,H,m,n);%Բ�壺4���߽�����+4���������壬mΪ��������ӳ�侫�ȵĲ���
%     ElementData=CirPlate2(n);%Բ�壺4����������
%    ElementData=ThreePyramid(a,b,h,n);%�壺3������
%     ElementData=ThreePyramidRevised(a,b,h,n);%�壺3������
%     ElementData=HybridPlane2(a,b,h,n);%�壺2������+2������
%     ElementData=HybridPlane(a,b,h,n);%�壺2������+2������+2������
    % ElementData=FiveTetra(a,b,h,n);%�壺5������
    % ElementData=DoubleTetraPyramid(a,b,h,n);%��������2������
    % ElementData=Hex_Prism(a,b,h,n);%�壺1������+2������
    % ElementData=DoubleHex(a,b,h,n);%�壺2������
    % ElementData=SingleHex(a,b,h,n);%�壺1������
%     ElementData=SinglePyramid(a,b,h,n);%��������1������
    % ElementData=DoublePrism(a,b,h,n);%�壺2������
    % ElementData=SingleTetra(a,b,h,n);%�����壺1������
    % ElementData=DoubleTetraPyramid2(a,b,h,n);%��������2������
    %###########################################################
    %% ��Ԫ������Ϣ��ȡ����Ԫ���ɶȵ�ȫ��������
    [ElementData,GlobalData] =Topology(ElementData);%
    %% ��Ԫ��������ģ��
    for i=1:length(ElementData)
        TYPE=ElementData{i}.Type;
        switch TYPE
            case 1
                [Ki,Mi] =EleMatrixHex(ElementData{i},E,mum,rhom);
            case 2
                [Ki,Mi] =EleMatrixPrism(ElementData{i},E,mum,rhom);
            case 3
                [Ki,Mi] =EleMatrixTetra(ElementData{i},E,mum,rhom);
            case 4
                [Ki,Mi] =EleMatrixPyram(ElementData{i},E,mum,rhom);
        end
        ElementData{i}.K=Ki;
        ElementData{i}.M=Mi;
    end
    %% ��Ԫ������װģ��
    [GK,GM]=Assemble(ElementData);
    GlobalData.GM=GM;
    GlobalData.GK=GK;
    %% �߽�����ģ��
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % �壺�ı߹�֧
%         FaceBC1={[1,0,0,0],[0,1,0,0],[1,0,0,-a],[0,1,0,-b]};%ʩ�ӱ߽�������λ��
%         bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC1,[1,1,1]);%�߽�������Ӧ�����ɶ�����
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % �壺�ı߼�֧
%         FaceBC1={[1,0,0,0],[1,0,0,-a]};
%         bc1=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC1,[0,1,1]);
%         FaceBC2={[0,1,0,0],[0,1,0,-b]};
%         bc2=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC2,[1,0,1]);
%         bc=bc1.*bc2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Բ�壺���ܹ�֧
%     FaceBCCP={[1,1,0,-100],[-1,1,0,-100],[-1,-1,0,-100],[1,-1,0,-100]}; bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBCCP,[1,1,1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %�������������֧
%     FaceBCCP={[0,0,1,0]};bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBCCP,[1,1,1]);
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % ����������֧
    FBC={[0,0,1,23]};bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FBC,[1,1,1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ����ģ�ͣ���ȫ����
%         bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH);solveway=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bn=find(bc);
    GK=GK(bn,bn);
    GM=GM(bn,bn);
    %%%%
    %% ����ֵ���ģ��
    [d,v,N]=SolveEigL(GK,GM,20,solveway);
    [sd,I]=sort(d);
    V=zeros(size(v));
    for i=1:length(I)
        V(:,i)=v(:,I(i));
    end
    Lmd=abs(sd).^0.5./(2*pi);
    D=E*h^3/(12*(1-mum^2));
    Lmdp=(sqrt(abs(sd))*b*b/(pi*pi))*sqrt(rhom*h/D);%����������Ƶ�ʲ���
    Lmdless=(2*rhom*(1+mum)*abs(sd)/E).^0.5*R;%Բ��������Ƶ�ʲ���
%% ���д��excel,��ֹ���׼���ʱ���Ա�������ʧ�����Ѿ�������ͽ׵�����
    Result=Lmd(1:12)';
%     Result=Lmd(1:12)';
%     Result=Lmdless(1:12)';
    str=['A',num2str(n)];
    xlswrite('result_Test.xlsx',Result,'Sheet1',str)
%%    
    %% �������ӻ�ģ��
    % ���񿪹�
    MeshSwitch=0;
    % �ڵ㿪��
    NodeSwitch=0;
    % �߽���������
    BCSwitch=0;
    % ģ̬����
    ModeSwitch=0;
    % λ�ƿ���
    
    if MeshSwitch==1
        VisualizeMesh(ElementData);%��Ԫ��ɫ�ɸ���
    end
    if NodeSwitch==1
        VisualizeNode(ElementData);
    end
    if BCSwitch==1
        VisualizeBC(ElementData,GlobalData.GlobalCoord,bc);
    end
    
    if ModeSwitch==1
        coef=50;%ģ̬λ�ƷŴ�ϵ��
        for i=1:7%�״�
            u=coef*V(:,i);
            ElementData=VectorRecovery(u,bc,ElementData);
            VisualizeMode(30,-1,ElementData);
        end
    end
end
