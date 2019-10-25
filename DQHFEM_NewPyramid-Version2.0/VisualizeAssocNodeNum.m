function VisualizeAssocNodeNum(ElementData,AssocNum)
%% 用于可视化公共节点的位置，便于调试单元组装是否出错


%% 单元节点可视化
Nele=length(ElementData);
%% 单元几何位置可视化
figure;
n=10;
for i=1:Nele
    MR=ElementData{i}.MapRight;
    Type=ElementData{i}.Type;
    switch Type
        case 1
            [r,s,t,tri]=StdHexPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3},tri{n+4}};
            if length(MR{1})==8
                [X,Y,Z]=MapCurHex(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurHex(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            C=ones(length(X),1);
        case 2
            [r,s,t,tri]=StdPrismPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
            if length(MR{1})==6
                [X,Y,Z]=MapCurPrism(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPrism(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
        case 3
            [r,s,t,tri]=StdTetraPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2}};
            if length(MR{1})==4
                [X,Y,Z]=MapCurTetra(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurTetra(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
        case 4
            [r,s,t,tri]=StdPyramPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
            if length(MR{1})==5
                [X,Y,Z]=MapCurPyram(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPyram(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
    end
    C=i*ones(length(X),1);
    for j=1:length(tri)
        trisurf(tri{j},X,Y,Z,C);
        alpha(0.2)
        shading interp
        hold on
    end
end


for i=1:length(AssocNum)
    nn=AssocNum{i};
    scatter3(ElementData{nn(1)}.X(nn(2)),ElementData{nn(1)}.Y(nn(2)),ElementData{nn(1)}.Z(nn(2)),20,'filled','r');
end

axis equal
axis off
view(3);
end
% Nele=length(ElementData);
% %% 全局节点坐标编号可视化
% figure
% for n=1:Nele
%     Vertex=ElementData{n}.MapRight;
%     Type=ElementData{n}.Type;
%     switch Type
%         case 1
%             PatchHex(Vertex);%单元可视化
%             hold on
%         case 2
%             PatchPrism(Vertex);%单元可视化
%             hold on
%         case 3
%             PatchTetra(Vertex);%单元可视化
%         case 4
%             PatchPyram(Vertex);%单元可视化
%             hold on
%     end
% end
% hold on
% GX=GlobalCoord{1};GY=GlobalCoord{2};GZ=GlobalCoord{3};
% GIH=[];
% for ni=1:Nele
%     GIH=[GIH,ElementData{ni}.GlobIndexH];
% end
% t=1;
% for i=1:length(GX)
%     if any(GIH==i)
%         t=t+1;
%     else
%         scatter3(GX(i),GY(i),GZ(i),20,'filled','b');
%         
%         text(GX(t),GY(t),GZ(t),num2str(t));
%         t=t+1;
%     end
% end
% axis equal
% view(3);
% end