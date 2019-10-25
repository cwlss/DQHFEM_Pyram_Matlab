function [NewElementData,GlobalData] =Topology(ElementData)
% 计算单元节点之间的拓扑关系（是否相邻节点）方便单元组装
Nele=length(ElementData);%读入单元个数
%% 计算表面节点坐标X,Y,Z以及表面节点个数、基函数个数；并将之存入相应单元结构体ElementData中
for n=1:Nele
    number=ElementData{n}.Number;
    MR=ElementData{n}.MapRight;
    Type=ElementData{n}.Type;
    switch Type
        case 1
            [Nbase,Nnode]=CountBasenNum(number,Type);
            [R,S,T] =FeketeHex(number);
            [XX,YY,ZZ]=MapCurHex(R,S,T,MR{1},MR{2},MR{3});%化曲为直，便于组装，节点判定时不会因为插值误差而出现公共节点判定错误
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
    %% 计算总体编号对应坐标GlobCood、单元索引号GlobIndex、单元内部索引号GlobIndexH
    %1.初始化每个单元的全部索引号GlobIndex/仅内部索引号GlobIndexH（每个单元均初始化为局部编号，首个单元准确，后续单元仍需判断修改）
    GlobIndex=cell(Nele,1);%边界+内部 自由度索引
    GlobIndexH=cell(Nele,1);%仅内部 自由度索引
    for n=1:Nele
        GlobIndex{n}=1:ElementData{n}.Nbase;
        GlobIndexH{n}=ElementData{n}.Nnode+1:ElementData{n}.Nbase;
    end
    ElementData{1}.GlobIndex=GlobIndex{1};
    ElementData{1}.GlobIndexH=GlobIndexH{1};
    %2.初始化全局编号的坐标值（长度待定，先将首个单元的节点坐标值作为全局节点坐标初始值）
    %因此，将第一个单元的节点坐标赋值给全局节点坐标
    Nh=ElementData{1}.Nbase-ElementData{1}.Nnode;
    GX=ElementData{1}.X;GX=[GX;1e15*ones(Nh,1)];%内部自由度不用施加边界条件，所以施加一个数值很大的虚拟坐标，以保证施加边界判断程序不会判断到该点
    GY=ElementData{1}.Y;GY=[GY;1e15*ones(Nh,1)];
    GZ=ElementData{1}.Z;GZ=[GZ;1e15*ones(Nh,1)];
    GXX=ElementData{1}.XX;GXX=[GXX;1e15*ones(Nh,1)];%内部自由度不用施加边界条件，所以施加一个数值很大的虚拟坐标，以保证施加边界判断程序不会判断到该点
    GYY=ElementData{1}.YY;GYY=[GYY;1e15*ones(Nh,1)];
    GZZ=ElementData{1}.ZZ;GZZ=[GZZ;1e15*ones(Nh,1)];
    
    
    p=0;%关联自由度计数器
    gn=ElementData{1}.Nbase;%独立自由度计数器(第1个单元默认独立。单元边界的自由度不一定独立，相互关联的多个自由度只算一个，内部的一定独立）
    for ni=2:Nele%遍历第ni个单元
        Np=ElementData{ni}.Nnode;%读入第ni个单元的边界自由度数量
        Nb=ElementData{ni}.Nbase;%读入第ni个单元的总自由度数量
        Nh=Nb-Np;%读入第ni个单元的内部（阶谱）自由度数量
        for i=1:Np
            %% 判断第ni个单元的第i个自由度是否与前面的所有单元的所有节点关联；
            % bool：1——共用
            % bool：0——独立
            bool=0;
            for nj=1:ni-1
                for j=1:ElementData{nj}.Nnode
                    r=(ElementData{ni}.XX(i)-ElementData{nj}.XX(j))^2+(ElementData{ni}.YY(i)-ElementData{nj}.YY(j))^2+(ElementData{ni}.ZZ(i)-ElementData{nj}.ZZ(j))^2;
                    if r<0.00000001
                        bool=1;
                        break%找到关联节点，跳出第一层循环
                    end
                end
                if bool
                    break%同时跳出第二层循环
                end
            end%判断结束
            %% 给第ni个单元的第i个元素赋予全局索引
            if bool==1
                p=p+1;%关联节点计数器+1
                GlobIndex{ni}(i)=GlobIndex{nj}(j);%关联节点：索引前后关联（后者关联前者）
                AssocNum{p}=[ni,i];
            else
                gn=gn+1;%独立节点计数器+1
                GlobIndex{ni}(i)=gn;%独立节点：建立新的全局索引编号
                %% 将独立节点坐标存入全局编号坐标GX,GY,GZ中
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
        GlobIndex{ni}(Np+1:Nb)=gn+1:gn+Nh;%补全GlobIndex中的内部索引
        GlobIndexH{ni}=gn+1:gn+Nh;
        ElementData{ni}.GlobIndex=GlobIndex{ni};%将单元索引GlobIndex集成到ElementData中
        ElementData{ni}.GlobIndexH=GlobIndexH{ni};%将单元内部索引GlobIndexH集成到ElementData中
        gn=gn+Nh;%建立内部自由度编号
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