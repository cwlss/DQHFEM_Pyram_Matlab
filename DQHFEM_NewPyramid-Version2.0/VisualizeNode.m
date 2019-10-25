function VisualizeNode(ElementData)
%% 模型所有的节点可视化



Nele=length(ElementData);
L=zeros(Nele,1);
for i=1:Nele
    L(i)=length(ElementData{i}.MapRight{1});
end
VisType=max(L);
%% 单元几何位置可视化
figure;
n=30;
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
            Vertex={MR{1}(1:8),MR{2}(1:8),MR{3}(1:8)};
        case 2
            [r,s,t,tri]=StdPrismPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
            if length(MR{1})==6
                [X,Y,Z]=MapCurPrism(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPrism(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            Vertex={MR{1}(1:6),MR{2}(1:6),MR{3}(1:6)};
        case 3
            [r,s,t,tri]=StdTetraPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2}};
            if length(MR{1})==4
                [X,Y,Z]=MapCurTetra(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurTetra(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            Vertex={MR{1}(1:4),MR{2}(1:4),MR{3}(1:4)};
        case 4
            [r,s,t,tri]=StdPyramPoint(n);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
            if length(MR{1})==5
                [X,Y,Z]=MapCurPyram(r,s,t,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPyram(r,s,t,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            Vertex={MR{1}(1:5),MR{2}(1:5),MR{3}(1:5)};
    end
if VisType<9
    PatchMesh(Vertex);
    hold on
else
    C=i*ones(length(X),1);
    for j=1:length(tri)
        trisurf(tri{j},X,Y,Z,C);
        alpha(0.2)
        shading interp 
        hold on
    end
end
for k=1:length(ElementData)
    scatter3(ElementData{k}.X,ElementData{k}.Y,ElementData{k}.Z,40,'filled','r');%节点可视化
    hold on
end
axis equal
axis off
view(3);
end