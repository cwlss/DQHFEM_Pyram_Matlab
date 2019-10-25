function PatchMesh(Vertex)
% 直边单元模型的可视化
N=length(Vertex{1});
switch N
    case 8
        C=[0,0,1];
        for i=1:6
            switch i
                case 1
                    pe=[1,2,3,4];
                case 2
                    pe=[1,2,6,5];
                case 3
                    pe=[4,3,7,8];
                case 4
                    pe=[4,1,5,8];
                case 5
                    pe=[5,6,7,8];
                case 6
                    pe=[2,3,7,6];
            end
            
            X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3));Vertex{1}(pe(4))];
            Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3));Vertex{2}(pe(4))];
            Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3));Vertex{3}(pe(4))];
            patch(X,Y,Z,C,'LineWidth',1.5);
            alpha(0.2);
        end
        axis off
    case 6
        Cq=[0.62745,0.12549,0.94118];
        Ct=[0.62745,0.12549,0.94118];
        for i=1:5
            switch i
                case 1
                    pe=[1,2,3];
                    X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3))];
                    Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3))];
                    Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3))];
                    patch(X,Y,Z,Ct,'LineWidth',1.5);
                    alpha(0.2);
                case 2
                    pe=[1,2,5,4];
                    X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3));Vertex{1}(pe(4))];
                    Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3));Vertex{2}(pe(4))];
                    Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3));Vertex{3}(pe(4))];
                    patch(X,Y,Z,Cq,'LineWidth',1.5);
                    alpha(0.2);
                case 3
                    pe=[2,3,6,5];
                    X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3));Vertex{1}(pe(4))];
                    Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3));Vertex{2}(pe(4))];
                    Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3));Vertex{3}(pe(4))];
                    patch(X,Y,Z,Cq,'LineWidth',1.5);
                    alpha(0.2);
                case 4
                    pe=[3,1,4,6];
                    X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3));Vertex{1}(pe(4))];
                    Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3));Vertex{2}(pe(4))];
                    Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3));Vertex{3}(pe(4))];
                    patch(X,Y,Z,Cq,'LineWidth',1.5);
                    alpha(0.2);
                case 5
                    pe=[4,5,6];
                    X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3))];
                    Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3))];
                    Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3))];
                    patch(X,Y,Z,Ct,'LineWidth',1.5);
                    alpha(0.2);
            end
        end
        axis off
    case 4
        C=[0,1,0];
        for k=1:4
            switch k
                case 1
                    pe=[1,2,3];
                case 2
                    pe=[1,2,4];
                case 3
                    pe=[2,3,4];
                case 4
                    pe=[3,1,4];
            end
            X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3))];
            Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3))];
            Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3))];
            patch(X,Y,Z,C,'LineWidth',1.5);
            alpha(0.2);
        end
        axis off
    case 5
        Ct=[1,1,1];
        Cq=[1,1,1];
        pe=[1,2,3,4];
        X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3));Vertex{1}(pe(4))];
        Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3));Vertex{2}(pe(4))];
        Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3));Vertex{3}(pe(4))];
        patch(X,Y,Z,Cq,'LineWidth',1.5);
        alpha(0.5);
        
        for i=2:5
            switch i
                case 2
                    pe=[1,2,5];
                case 3
                    pe=[2,3,5];
                case 4
                    pe=[3,4,5];
                case 5
                    pe=[4,1,5];
            end
            
            X=[Vertex{1}(pe(1));Vertex{1}(pe(2));Vertex{1}(pe(3))];
            Y=[Vertex{2}(pe(1));Vertex{2}(pe(2));Vertex{2}(pe(3))];
            Z=[Vertex{3}(pe(1));Vertex{3}(pe(2));Vertex{3}(pe(3))];
            patch(X,Y,Z,Ct,'LineWidth',1.5);
            alpha(0.5);
        end
        axis off
        axis equal
end


