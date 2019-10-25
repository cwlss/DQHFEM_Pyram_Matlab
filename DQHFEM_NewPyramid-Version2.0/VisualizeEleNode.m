function VisualizeEleNode(ElementData,GlobalCoord)
%% 每个单元上的节点，分步可视化


%% 单个单元全局坐标可视化
Nele=length(ElementData);
GX=GlobalCoord{1};GY=GlobalCoord{2};GZ=GlobalCoord{3};
for i=1:Nele
    figure;
    axis equal
    for n=1:Nele
        Vertex=ElementData{n}.MapRight;
        Type=ElementData{n}.Type;
        switch Type
            case 1
                PatchHex(Vertex);%单元可视化
                hold on
            case 2
                PatchPrism(Vertex);%单元可视化
                hold on
            case 3
                PatchTetra(Vertex);%单元可视化
            case 4
                PatchPyram(Vertex);%单元可视化
                hold on
        end
    end
    view(3)
    hold on
    for ii=1:ElementData{i}.Nnode
        num=ElementData{i}.GlobIndex(ii);
        scatter3(GX(num),GY(num),GZ(num),'filled','r');
        text(GX(num),GY(num),GZ(num),num2str(num));
    end
end
end