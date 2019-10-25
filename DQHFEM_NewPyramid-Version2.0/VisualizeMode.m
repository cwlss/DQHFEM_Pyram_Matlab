function ElementData=VisualizeMode(n,direction,ElementData)
%% 模态可视化
% n 为场图像精细度，每个方向布点个数，三个方向相同
% direction：1:U 2:V 3:W 0:无变形：-1：真实位移
figure
Nele=length(ElementData);
for i=1:Nele
    type=ElementData{i}.Type;
    MR=ElementData{i}.MapRight;
    number=ElementData{i}.Number;
    U=ElementData{i}.eleU;
    V=ElementData{i}.eleV;
    W=ElementData{i}.eleW;
    n=n(1);
    switch type
        case 1
            [R,S,T,tri]=StdHexPoint(n);
            if length(MR{1})==8
                [X,Y,Z]=MapCurHex(R,S,T,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurHex(R,S,T,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            [G,Gr,Gs,Gt]=BaseHex(R,S,T,number);
            T=InterpHex(number);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3},tri{n+4}};
        case 2
            [R,S,T,tri]=StdPrismPoint(n);
            if length(MR{1})==6
                [X,Y,Z]=MapCurPrism(R,S,T,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPrism(R,S,T,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            [G,Gr,Gs,Gt]=BasePrism(R,S,T,number);
            T=InterpPrism(number);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
        case 3
            [R,S,T,tri]=StdTetraPoint(n);
            if length(MR{1})==4
                [X,Y,Z]=MapCurTetra(R,S,T,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurTetra(R,S,T,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            [G,Gr,Gs,Gt]=BaseTetra(R,S,T,number);
            T=InterpTetra(number);
            tri={tri{1},tri{n},tri{n+1},tri{n+2}};
        case 4
            [R,S,T,tri]=StdPyramPoint(n);
            if length(MR{1})==5
                [X,Y,Z]=MapCurPyram(R,S,T,MR{1},MR{2},MR{3});
            else
                [X,Y,Z]=MapCurPyram(R,S,T,MR{1},MR{2},MR{3},ElementData{i}.MapNumber);
            end
            [G,Gr,Gs,Gt]=BasePyram(R,S,T,number);
            T=InterpPyram(number);
            tri={tri{1},tri{n},tri{n+1},tri{n+2},tri{n+3}};
    end
    [GG,~,~,~]=InterpTrans(G,Gr,Gs,Gt,T);
    u=GG*U;v=GG*V;w=GG*W;
    uvw=sqrt(u.^2+v.^2+w.^2);
    switch direction
        case 1
            C=u;
        case 2
            C=v;
        case 3
            C=w;
        case 0
            C=zeros(size(u));
        case -1
            C=uvw;
    end
    x=X+u;
    y=Y+v;
    z=Z+w;
    for j=1:length(tri)
        trisurf(tri{j},x,y,z,C)
        alpha(1)
        hold on
    end
    hold on
    axis equal;
    axis off
    shading interp
    %     light
    %     lighting flat
    %     lighting gouraud
    %     lighting phong
    %     lighting none
    ElementData{i}.XField=X;
    ElementData{i}.YField=Y;
    ElementData{i}.ZField=Z;
    ElementData{i}.UField=u;
    ElementData{i}.VField=v;
    ElementData{i}.WField=w;
end
end



