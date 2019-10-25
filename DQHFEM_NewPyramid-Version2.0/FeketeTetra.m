function [R,S,T] =FeketeTetra(number)

%% 顶点
R=[0;1;0;0];
S=[0;0;1;0];
T=[0;0;0;1];
%% 边内结点
for i=1:6
    ne=number.edge(i);
    x=GaussLobattoR(ne+2,0,1);
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
    x=GaussLobattoR(nf+3,0,1);
    [X,Y]=TrigLobatto(x);
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
function [S, T, Se, Te]=TrigLobatto(t)
%
% TrigLobatto: Generate non-equally spaces nodes on a unit right triangle
%
% Calling Sequence:
%
%     [S, T]=TrigLobatto(t)
%
%     [S, T, Se, Te]=TrigLobatto(t)
%
% INPUT:
%
%     t :  An unequally space node vector defined on [0, 1]
%
% OUTPUT:
%
%    S, T :  Natrual coordinates of unequally spaced nodes
%            on a unit right triangle.
%
%    Se, Te :  Natrual coordinates of equally spaced nodes
%            on a unit right triangle.
%

m=length(t);
s=linspace(0, 1, m);
N=m*(m+1)/2;

% Equal nodes on unit right triangle
if nargout>2
    Se=zeros(N, 1); Te=Se;
end
Is=zeros(N, 1); It=Is;
p=1;
for j=1:m
    for i=1:(m-j)+1
        if nargout>2
            Se(p)=s(i); Te(p)=s(j);
        end
        Is(p)=i; It(p)=j;
        p=p+1;
    end
end

% Transform to be non-equally spaced
Ir=m+2-Is-It;
Sn=t(Is); Tn=t(It); Rn=t(Ir);
Pn=9*Sn.*Tn.*Rn;
Dn=Sn+Tn+Rn+Pn;
S=(Sn+Pn/3)./Dn; T=(Tn+Pn/3)./Dn;
end

% test successful!!!




