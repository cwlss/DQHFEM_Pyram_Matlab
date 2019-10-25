function [R,S,T] =UniformPrism(number)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
%% 顶点
R=[0;1;0;0;1;0];
S=[0;0;1;0;0;1];
T=[0;0;0;1;1;1];
%% 边内结点
for i=1:9
    ne=number.edge(i);
    x=linspace(0,1,ne+2);x=x';
    switch i
        case 1
            Re=x(2:ne+1);Se=zeros(ne,1);Te=Se;
        case 2
            Re=x(2:ne+1);Se=1-x(2:ne+1);Te=zeros(ne,1);
        case 3
            Re=zeros(ne,1);Se=x(2:ne+1);Te=zeros(ne,1);
        case 4
            Re=zeros(ne,1);Se=zeros(ne,1);Te=x(2:ne+1);
        case 5
            Re=ones(ne,1);Se=zeros(ne,1);Te=x(2:ne+1);
        case 6
            Re=zeros(ne,1);Se=ones(ne,1);Te=x(2:ne+1);
        case 7
            Re=x(2:ne+1);Se=zeros(ne,1);Te=ones(ne,1);
        case 8
            Re=x(2:ne+1);Se=1-x(2:ne+1);Te=ones(ne,1);
        case 9
            Re=zeros(ne,1);Se=x(2:ne+1);Te=ones(ne,1);
    end
    R=[R;Re];S=[S;Se];T=[T;Te];
end
%% 面内节点
%% triface1
fi=1;
nf=number.face{fi};
Rf=zeros((nf+1)*nf/2,1);Sf=Rf;Tf=Rf;
x=linspace(0,1,nf+3);x=x';
[r, s]=TrigUniform(x);
t=zeros(size(r));
p=1;
for j=1:nf
    for i=1:nf-j+1
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf-j+4)/2+i+1;
        Rf(p)=r(nt);Sf(p)=s(nt);Tf(p)=t(nt);
        p=p+1;
    end
end
R=[R;Rf];S=[S;Sf];T=[T;Tf];
%% qface 2 3 4
for fi=2:4
    nf=number.face{fi};
    x1=linspace(0,1,nf(1)+2);x1=x1';x=x1(2:nf(1)+1);
    x2=linspace(0,1,nf(2)+2);x2=x2';y=x2(2:nf(2)+1);
    Rf=zeros(nf(1)*nf(2),1);Sf=Rf;Tf=Rf;
    switch fi
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
                    Rf(ij,1)=x(i);
                    Sf(ij,1)=1-x(i);
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
        case 4
            ij=1;
            for j=1:nf(2)
                for i=1:nf(1)
                    Rf(ij,1)=0;
                    Sf(ij,1)=x(i);
                    Tf(ij,1)=y(j);
                    ij=ij+1;
                end
            end
    end
    R=[R;Rf];S=[S;Sf];T=[T;Tf];
end
%% triface5
fi=5;
nf=number.face{fi};
Rf=zeros((nf+1)*nf/2,1);Sf=Rf;Tf=Rf;
x=linspace(0,1,nf+3);x=x';
[r, s]=TrigUniform(x);
t=ones(size(R));
p=1;
for j=1:nf
    for i=1:nf-j+1
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf-j+4)/2+i+1;
        Rf(p)=r(nt);Sf(p)=s(nt);Tf(p)=t(nt);
        p=p+1;
    end
end
R=[R;Rf];S=[S;Sf];T=[T;Tf];
end
function [S, T]=TrigUniform(t)

n=length(t);
N=n*(n+1)/2;
S=zeros(N,1);T=S;
ij=1;
for j=1:n
    for i=1:n-j+1
        S(ij)=t(i);
        T(ij)=t(j);
        ij=ij+1;
    end
end
end

% test successful!!!




