function [R,S,T] =UniformTetra(number)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
%% 顶点
R=[0;1;0;0];
S=[0;0;1;0];
T=[0;0;0;1];
%% 边内结点
for i=1:6
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
            Re=x(2:ne+1);Se=zeros(ne,1);Te=1-x(2:ne+1);
        case 6
            Re=zeros(ne,1);Se=x(2:ne+1);Te=1-x(2:ne+1);
    end
    R=[R;Re];S=[S;Se];T=[T;Te];
end
%% 面内节点
%% triface 1 2 3 4
for fi=1:4;
    nf=number.face{fi};
    Rf=zeros((nf+1)*nf/2,1);Sf=Rf;Tf=Rf;
    x=linspace(0,1,nf+3);x=x';
    [X,Y]=TrigUniform(x);
    switch fi
        case 1
            r=X;s=Y;t=zeros(size(r));
        case 2
            r=X;s=zeros(size(r));t=Y;
        case 3
            r=X;s=Y;t=1-X-Y;
        case 4
            r=zeros(size(r));s=X;t=Y;
            
    end
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




