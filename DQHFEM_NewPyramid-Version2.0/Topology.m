function [NewElementData,GlobalData] =Topology(ElementData)
% ���㵥Ԫ�ڵ�֮������˹�ϵ���Ƿ����ڽڵ㣩���㵥Ԫ��װ
Nele=length(ElementData);%���뵥Ԫ����
%% �������ڵ�����X,Y,Z�Լ�����ڵ����������������������֮������Ӧ��Ԫ�ṹ��ElementData��
for n=1:Nele
    number=ElementData{n}.Number;
    MR=ElementData{n}.MapRight;
    Type=ElementData{n}.Type;
    switch Type
        case 1
            [Nbase,Nnode]=CountBasenNum(number,Type);
            [R,S,T] =FeketeHex(number);
            [XX,YY,ZZ]=MapCurHex(R,S,T,MR{1},MR{2},MR{3});%����Ϊֱ��������װ���ڵ��ж�ʱ������Ϊ��ֵ�������ֹ����ڵ��ж�����
            if length(MR{1})==8
                X=XX;Y=YY;Z=ZZ;
            else
                [X,Y,Z]=MapCurHex(R,S,T,MR{1},MR{2},MR{3},ElementData{n}.MapNumber);
            end
        case 2
            [Nbase,Nnode]=CountBasenNum(number,Type);
            [R,S,T] =FeketePrism(number);
            [XX,YY,ZZ]=MapCurPrism(R,S,T,MR{1},MR{2},MR{3});
            if length(MR{1})==6
                X=XX;Y=YY;Z=ZZ;
            else
                [X,Y,Z]=MapCurPrism(R,S,T,MR{1},MR{2},MR{3},ElementData{n}.MapNumber);
            end
        case 3
            [Nbase,Nnode]=CountBasenNum(number,Type);
            [R,S,T] =FeketeTetra(number);
            [XX,YY,ZZ]=MapCurTetra(R,S,T,MR{1},MR{2},MR{3});
            if length(MR{1})==4
                X=XX;Y=YY;Z=ZZ;
            else
                [X,Y,Z]=MapCurTetra(R,S,T,MR{1},MR{2},MR{3},ElementData{n}.MapNumber);
            end
        case 4
            [Nbase,Nnode]=CountBasenNum(number,Type);
            [R,S,T] =FeketePyram(number);
            [XX,YY,ZZ]=MapCurPyram(R,S,T,MR{1},MR{2},MR{3});
            if length(MR{1})==5
                X=XX;Y=YY;Z=ZZ;
            else
                [X,Y,Z]=MapCurPyram(R,S,T,MR{1},MR{2},MR{3},ElementData{n}.MapNumber);
            end
    end
    ElementData{n}.X=X;ElementData{n}.Y=Y;ElementData{n}.Z=Z;
    ElementData{n}.XX=XX;ElementData{n}.YY=YY;ElementData{n}.ZZ=ZZ;
    ElementData{n}.Nbase=Nbase;ElementData{n}.Nnode=Nnode;
end

if length(ElementData)==1
    Xh=1e15*ones(ElementData{1}.Nbase-ElementData{1}.Nnode,1);
    ElementData{1}.GlobIndex=1:ElementData{1}.Nbase;
    ElementData{1}.GlobIndexH=ElementData{1}.Nnode+1:ElementData{1}.Nbase;
    GX=[ElementData{1}.X;Xh];
    GY=[ElementData{1}.Y;Xh];
    GZ=[ElementData{1}.Z;Xh];
    GXX=[ElementData{1}.XX;Xh];
    GYY=[ElementData{1}.YY;Xh];
    GZZ=[ElementData{1}.ZZ;Xh];
    GlobalCoord={GX,GY,GZ};
    LinGlobalCoord={GXX,GYY,GZZ};
    ElementData{1}.NG=length(GX);
    AssocNum=[];
    NewElementData=ElementData;
    GlobalData.GlobalCoord=GlobalCoord;
    GlobalData.LinGlobalCoord=LinGlobalCoord;
    GlobalData.AssocNum=AssocNum;
    GlobalData.IndexH=ElementData{1}.GlobIndexH;
else
    %% ���������Ŷ�Ӧ����GlobCood����Ԫ������GlobIndex����Ԫ�ڲ�������GlobIndexH
    %1.��ʼ��ÿ����Ԫ��ȫ��������GlobIndex/���ڲ�������GlobIndexH��ÿ����Ԫ����ʼ��Ϊ�ֲ���ţ��׸���Ԫ׼ȷ��������Ԫ�����ж��޸ģ�
    GlobIndex=cell(Nele,1);%�߽�+�ڲ� ���ɶ�����
    GlobIndexH=cell(Nele,1);%���ڲ� ���ɶ�����
    for n=1:Nele
        GlobIndex{n}=1:ElementData{n}.Nbase;
        GlobIndexH{n}=ElementData{n}.Nnode+1:ElementData{n}.Nbase;
    end
    ElementData{1}.GlobIndex=GlobIndex{1};
    ElementData{1}.GlobIndexH=GlobIndexH{1};
    %2.��ʼ��ȫ�ֱ�ŵ�����ֵ�����ȴ������Ƚ��׸���Ԫ�Ľڵ�����ֵ��Ϊȫ�ֽڵ������ʼֵ��
    %��ˣ�����һ����Ԫ�Ľڵ����긳ֵ��ȫ�ֽڵ�����
    Nh=ElementData{1}.Nbase-ElementData{1}.Nnode;
    GX=ElementData{1}.X;GX=[GX;1e15*ones(Nh,1)];%�ڲ����ɶȲ���ʩ�ӱ߽�����������ʩ��һ����ֵ�ܴ���������꣬�Ա�֤ʩ�ӱ߽��жϳ��򲻻��жϵ��õ�
    GY=ElementData{1}.Y;GY=[GY;1e15*ones(Nh,1)];
    GZ=ElementData{1}.Z;GZ=[GZ;1e15*ones(Nh,1)];
    GXX=ElementData{1}.XX;GXX=[GXX;1e15*ones(Nh,1)];%�ڲ����ɶȲ���ʩ�ӱ߽�����������ʩ��һ����ֵ�ܴ���������꣬�Ա�֤ʩ�ӱ߽��жϳ��򲻻��жϵ��õ�
    GYY=ElementData{1}.YY;GYY=[GYY;1e15*ones(Nh,1)];
    GZZ=ElementData{1}.ZZ;GZZ=[GZZ;1e15*ones(Nh,1)];
    
    
    p=0;%�������ɶȼ�����
    gn=ElementData{1}.Nbase;%�������ɶȼ�����(��1����ԪĬ�϶�������Ԫ�߽�����ɶȲ�һ���������໥�����Ķ�����ɶ�ֻ��һ�����ڲ���һ��������
    for ni=2:Nele%������ni����Ԫ
        Np=ElementData{ni}.Nnode;%�����ni����Ԫ�ı߽����ɶ�����
        Nb=ElementData{ni}.Nbase;%�����ni����Ԫ�������ɶ�����
        Nh=Nb-Np;%�����ni����Ԫ���ڲ������ף����ɶ�����
        for i=1:Np
            %% �жϵ�ni����Ԫ�ĵ�i�����ɶ��Ƿ���ǰ������е�Ԫ�����нڵ������
            % bool��1��������
            % bool��0��������
            bool=0;
            for nj=1:ni-1
                for j=1:ElementData{nj}.Nnode
                    r=(ElementData{ni}.XX(i)-ElementData{nj}.XX(j))^2+(ElementData{ni}.YY(i)-ElementData{nj}.YY(j))^2+(ElementData{ni}.ZZ(i)-ElementData{nj}.ZZ(j))^2;
                    if r<0.00000001
                        bool=1;
                        break%�ҵ������ڵ㣬������һ��ѭ��
                    end
                end
                if bool
                    break%ͬʱ�����ڶ���ѭ��
                end
            end%�жϽ���
            %% ����ni����Ԫ�ĵ�i��Ԫ�ظ���ȫ������
            if bool==1
                p=p+1;%�����ڵ������+1
                GlobIndex{ni}(i)=GlobIndex{nj}(j);%�����ڵ㣺����ǰ����������߹���ǰ�ߣ�
                AssocNum{p}=[ni,i];
            else
                gn=gn+1;%�����ڵ������+1
                GlobIndex{ni}(i)=gn;%�����ڵ㣺�����µ�ȫ���������
                %% �������ڵ��������ȫ�ֱ������GX,GY,GZ��
                GX=[GX;ElementData{ni}.X(i)];
                GY=[GY;ElementData{ni}.Y(i)];
                GZ=[GZ;ElementData{ni}.Z(i)];
                GXX=[GXX;ElementData{ni}.XX(i)];
                GYY=[GYY;ElementData{ni}.YY(i)];
                GZZ=[GZZ;ElementData{ni}.ZZ(i)];
            end
        end
        GX=[GX;1e15*ones(Nh,1)];
        GY=[GY;1e15*ones(Nh,1)];
        GZ=[GZ;1e15*ones(Nh,1)];
        GXX=[GXX;1e15*ones(Nh,1)];
        GYY=[GYY;1e15*ones(Nh,1)];
        GZZ=[GZZ;1e15*ones(Nh,1)];
        GlobIndex{ni}(Np+1:Nb)=gn+1:gn+Nh;%��ȫGlobIndex�е��ڲ�����
        GlobIndexH{ni}=gn+1:gn+Nh;
        ElementData{ni}.GlobIndex=GlobIndex{ni};%����Ԫ����GlobIndex���ɵ�ElementData��
        ElementData{ni}.GlobIndexH=GlobIndexH{ni};%����Ԫ�ڲ�����GlobIndexH���ɵ�ElementData��
        gn=gn+Nh;%�����ڲ����ɶȱ��
    end
    GlobalCoord={GX,GY,GZ};
    LinGlobalCoord={GXX,GYY,GZZ};
    ElementData{1}.NG=length(GX);
    NewElementData=ElementData;
    GlobalData.GlobalCoord=GlobalCoord;
    GlobalData.LinGlobalCoord=LinGlobalCoord;
    GlobalData.AssocNum=AssocNum;
    
    GIH=[];
    for ni=1:Nele
        GIH=[GIH,ElementData{ni}.GlobIndexH];
    end
    GlobalData.IndexH=GIH;
end
end