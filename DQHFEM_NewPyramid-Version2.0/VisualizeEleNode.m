function VisualizeEleNode(ElementData,GlobalCoord)
%% ÿ����Ԫ�ϵĽڵ㣬�ֲ����ӻ�


%% ������Ԫȫ��������ӻ�
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
                PatchHex(Vertex);%��Ԫ���ӻ�
                hold on
            case 2
                PatchPrism(Vertex);%��Ԫ���ӻ�
                hold on
            case 3
                PatchTetra(Vertex);%��Ԫ���ӻ�
            case 4
                PatchPyram(Vertex);%��Ԫ���ӻ�
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