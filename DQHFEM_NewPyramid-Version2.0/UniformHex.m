function [R,S,T] =UniformHex(number)
%% 与Fekete对立的均匀插值节点计算程序，可供选择

%% 顶点
R=[0;1;1;0;0;1;1;0];
S=[0;0;1;1;0;0;1;1];
T=[0;0;0;0;1;1;1;1];
%% 边内结点
for i=1:12
    ne=number.edge(i);
    %x=GaussLobattoR(ne+2,0,1);
    x=linspace(0,1,ne+2);x=x';
    switch i
        case 1
            Re=x(2:ne+1);Se=zeros(ne,1);Te=Se;
        case 2
            Re=ones(ne,1);Se=x(2:ne+1);Te=zeros(ne,1);
        case 3
            Re=x(2:ne+1);Se=ones(ne,1);Te=zeros(ne,1);
        case 4
            Re=zeros(ne,1);Se=x(2:ne+1);Te=zeros(ne,1);
        case 5
            Re=zeros(ne,1);Se=zeros(ne,1);Te=x(2:ne+1);
        case 6
            Re=ones(ne,1);Se=zeros(ne,1);Te=x(2:ne+1);
        case 7
            Re=ones(ne,1);Se=ones(ne,1);Te=x(2:ne+1);
        case 8
            Re=zeros(ne,1);Se=ones(ne,1);Te=x(2:ne+1);
        case 9
            Re=x(2:ne+1);Se=zeros(ne,1);Te=ones(ne,1);
        case 10
            Re=ones(ne,1);Se=x(2:ne+1);Te=ones(ne,1);
        case 11
            Re=x(2:ne+1);Se=ones(ne,1);Te=ones(ne,1);
        case 12
            Re=zeros(ne,1);Se=x(2:ne+1);Te=ones(ne,1);
    end
    R=[R;Re];S=[S;Se];T=[T;Te];
end
%% 面内节点
for fi=1:6
    nf=number.face{fi};
    x1=linspace(0,1,nf(1)+2);x=x1(2:nf(1)+1);x1=x1';
    x2=linspace(0,1,nf(2)+2);y=x2(2:nf(2)+1);x2=x2';
    Rf=zeros(nf(1)*nf(2),1);Sf=Rf;Tf=Rf;
    switch fi
        case 1
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=x(i);
                    Sf(ij,1)=y(j);
                    Tf(ij,1)=0;
                    ij=ij+1;
                end
            end
        case 2
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=x(i);
                    Sf(ij,1)=0;
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
        case 3
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=1;
                    Sf(ij,1)=x(i);
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
        case 4
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=x(i);
                    Sf(ij,1)=1;
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
        case 5
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=0;
                    Sf(ij,1)=x(i);
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
        case 6
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=x(i);
                    Sf(ij,1)=y(j);
                    Tf(ij,1)=1;
                    ij=ij+1;
                end
            end
    end
    R=[R;Rf];S=[S;Sf];T=[T;Tf];
end
end
% test successful!




