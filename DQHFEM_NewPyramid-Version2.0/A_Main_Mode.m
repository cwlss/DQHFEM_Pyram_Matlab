%%2019.07.23 by Guomao
clear;clc
%% 分析类型：模态
% 特征值求解算法
solveway=2;%部分求解法，不适用于自由边界条件
%% 前处理模块
% 材料参数
E=195E9;
mum=0.3;
rhom=7722.3;
% rhom=100;
% 基本尺寸
a=100;b=100;h=100;%方板
R=100;H=20;m=15;%圆板：m为曲边映射节点的单方向个数，m建议取20以下，太大会导致插值转换矩阵接近奇异（由于矩阵元素均小于1，行列式维度越大，值接近于0）
% 单元阶次
for n=2
    %##################################################################
    % ElementData元胞矩阵：长度为单元数量，元素为各个单元对应结构体，每个结构体共有的字段如下
    % 'Type'字段：单元类型信息
    % 'Number'字段：单元边界的基函数插值节点数量，与阶次n相关
    % 'MapRight'字段:映射权系数，为单元映射节点的几何坐标（注意有顺序要求）
    % 'MapNumber'字段：单元边界的几何映射节点数量，与边界选取的映射节点数目和位置相关
    
% 以下均为模型录入程序（前处理过程），取消注释即可更改算例进行计算，
% 如果需要新建立模型，重新建立一个录入程序即可，
% 以下录入程序若有不需要的，可从程序包中删掉，不会影响整个流程的计算

    ElementData=huojian(n);%火箭
% ElementData=CirclePlate1(R,H,m,n);%圆板映射不完全算例
%     ElementData=CirPlate1(R,H,m,n);%圆板：4曲边金字塔+4曲边四面体，m为控制曲边映射精度的参数
%     ElementData=CirPlate2(n);%圆板：4曲边三棱柱
%    ElementData=ThreePyramid(a,b,h,n);%板：3金字塔
%     ElementData=ThreePyramidRevised(a,b,h,n);%板：3金字塔
%     ElementData=HybridPlane2(a,b,h,n);%板：2金字塔+2四面体
%     ElementData=HybridPlane(a,b,h,n);%板：2三棱柱+2金字塔+2四面体
    % ElementData=FiveTetra(a,b,h,n);%板：5四面体
    % ElementData=DoubleTetraPyramid(a,b,h,n);%金字塔：2四面体
    % ElementData=Hex_Prism(a,b,h,n);%板：1六面体+2三棱柱
    % ElementData=DoubleHex(a,b,h,n);%板：2六面体
    % ElementData=SingleHex(a,b,h,n);%板：1六面体
%     ElementData=SinglePyramid(a,b,h,n);%金字塔：1金字塔
    % ElementData=DoublePrism(a,b,h,n);%板：2三棱柱
    % ElementData=SingleTetra(a,b,h,n);%四面体：1四面体
    % ElementData=DoubleTetraPyramid2(a,b,h,n);%金字塔：2四面体
    %###########################################################
    %% 单元拓扑信息读取（单元自由度的全局索引）
    [ElementData,GlobalData] =Topology(ElementData);%
    %% 单元矩阵生成模块
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
    %% 单元矩阵组装模块
    [GK,GM]=Assemble(ElementData);
    GlobalData.GM=GM;
    GlobalData.GK=GK;
    %% 边界条件模块
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 板：四边固支
%         FaceBC1={[1,0,0,0],[0,1,0,0],[1,0,0,-a],[0,1,0,-b]};%施加边界条件的位置
%         bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC1,[1,1,1]);%边界条件对应的自由度索引
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % 板：四边简支
%         FaceBC1={[1,0,0,0],[1,0,0,-a]};
%         bc1=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC1,[0,1,1]);
%         FaceBC2={[0,1,0,0],[0,1,0,-b]};
%         bc2=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBC2,[1,0,1]);
%         bc=bc1.*bc2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %圆板：四周固支
%     FaceBCCP={[1,1,0,-100],[-1,1,0,-100],[-1,-1,0,-100],[1,-1,0,-100]}; bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBCCP,[1,1,1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %金字塔：底面固支
%     FaceBCCP={[0,0,1,0]};bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FaceBCCP,[1,1,1]);
%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 火箭：底面固支
    FBC={[0,0,1,23]};bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH,FBC,[1,1,1]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 任意模型：完全自由
%         bc=ApplyFaceBC(GlobalData.LinGlobalCoord,GlobalData.IndexH);solveway=1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bn=find(bc);
    GK=GK(bn,bn);
    GM=GM(bn,bn);
    %%%%
    %% 特征值求解模块
    [d,v,N]=SolveEigL(GK,GM,20,solveway);
    [sd,I]=sort(d);
    V=zeros(size(v));
    for i=1:length(I)
        V(:,i)=v(:,I(i));
    end
    Lmd=abs(sd).^0.5./(2*pi);
    D=E*h^3/(12*(1-mum^2));
    Lmdp=(sqrt(abs(sd))*b*b/(pi*pi))*sqrt(rhom*h/D);%方板的无因次频率参数
    Lmdless=(2*rhom*(1+mum)*abs(sd)/E).^0.5*R;%圆板的无因次频率参数
%% 结果写入excel,防止升阶计算时电脑崩溃，遗失所有已经计算过低阶的数据
    Result=Lmd(1:12)';
%     Result=Lmd(1:12)';
%     Result=Lmdless(1:12)';
    str=['A',num2str(n)];
    xlswrite('result_Test.xlsx',Result,'Sheet1',str)
%%    
    %% 后处理：可视化模块
    % 网格开关
    MeshSwitch=0;
    % 节点开关
    NodeSwitch=0;
    % 边界条件开关
    BCSwitch=0;
    % 模态开关
    ModeSwitch=0;
    % 位移开关
    
    if MeshSwitch==1
        VisualizeMesh(ElementData);%单元颜色可更改
    end
    if NodeSwitch==1
        VisualizeNode(ElementData);
    end
    if BCSwitch==1
        VisualizeBC(ElementData,GlobalData.GlobalCoord,bc);
    end
    
    if ModeSwitch==1
        coef=50;%模态位移放大系数
        for i=1:7%阶次
            u=coef*V(:,i);
            ElementData=VectorRecovery(u,bc,ElementData);
            VisualizeMode(30,-1,ElementData);
        end
    end
end
