function   [R,S,T]=FeketePyram(number)
%%
%输入：number： 和基函数同步的Fekete点数量
%输出：R,S,T：  Fekete点在三个方向上的坐标列向量
%% qface
nqr=number.face{1}(1);
nqs=number.face{1}(2);
x=GaussLobattoR(nqr+2,-1,1);
y=GaussLobattoR(nqs+2,-1,1);
N=nqr*nqs;
Rfq=zeros(N,1);Sfq=Rfq;Tfq=Rfq;
p=1;
for j=1:nqs
    for i=1:nqr
        Rfq(p)=x(i+1);Sfq(p)=y(j+1);
        p=p+1;
    end
end
%% triface1
nf=number.face{2};
Rf1=zeros((nf+1)*nf/2,1);Sf1=Rf1;Tf1=Rf1;
x=GaussLobattoR(nf+3,0,1);
[r, s]=TrigLobatto(x);
% [r,s]=TrigNodeVect(nf+3);

R1=2*r+s-1;
S1=s-1;
T1=s;
p=1;
for j=1:nf
    for i=1:nf-j+1
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf-j+4)/2+i+1;
        Rf1(p)=R1(nt);Sf1(p)=S1(nt);Tf1(p)=T1(nt);
        p=p+1;
    end
end
%% triface2
nf=number.face{3};
Rf2=zeros((nf+1)*nf/2,1);Sf2=Rf2;Tf2=Rf2;
x=GaussLobattoR(nf+3,0,1);
[r, s]=TrigLobatto(x);
% [r,s]=TrigNodeVect(nf+3);

R2=1-s;
S2=2*r+s-1;
T2=s;
p=1;
for j=1:nf
    for i=1:nf-j+1
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf+3-j+1)/2+i+1;
        Rf2(p)=R2(nt);Sf2(p)=S2(nt);Tf2(p)=T2(nt);
        p=p+1;
    end
end
%% triface3
nf=number.face{4};
Rf3=zeros((nf+1)*nf/2,1);Sf3=Rf3;Tf3=Rf3;
x=GaussLobattoR(nf+3,0,1);
[r, s]=TrigLobatto(x);
% [r,s]=TrigNodeVect(nf+3);

R3=1-2*r-s;
S3=1-s;
T3=s;
p=1;
for j=1:nf
    for i=1:nf+1-j
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf-j+4)/2+i+1;
        Rf3(p)=R3(nt);Sf3(p)=S3(nt);Tf3(p)=T3(nt);
        p=p+1;
    end
end
%% triface4
nf=number.face{5};
Rf4=zeros((nf+1)*nf/2,1);Sf4=Rf4;Tf4=Rf4;
x=GaussLobattoR(nf+3,0,1);
[r, s]=TrigLobatto(x);
% [r,s]=TrigNodeVect(nf+3);

R4=s-1;
S4=1-2*r-s;
T4=s;
p=1;
for j=1:nf
    for i=1:nf-j+1
        nt=(nf+3)*(nf+4)/2-(nf+3-j)*(nf-j+4)/2+i+1;
        Rf4(p)=R4(nt);Sf4(p)=S4(nt);Tf4(p)=T4(nt);
        p=p+1;
    end
end
%% edge1
ne=number.edge(1);
x=GaussLobattoR(ne+2,0,1);
Re1=2*x(2:ne+1)-1;
Se1=-1*ones(size(Re1));
Te1=0*Re1;
%% edge2
ne=number.edge(2);
x=GaussLobattoR(ne+2,0,1);
Se2=2*x(2:ne+1)-1;
Re2=ones(size(Se2));
Te2=0*Se2;
%% edge3
ne=number.edge(3);
x=GaussLobattoR(ne+2,0,1);
Re3=2*x(ne+1:-1:2)-1;
Se3=1*ones(size(Re3));
Te3=0*Re3;
%% edge4
ne=number.edge(4);
x=GaussLobattoR(ne+2,0,1);
Se4=2*x(ne+1:-1:2)-1;
Re4=-1*ones(size(Se4));
Te4=0*Se4;
%% edge5
ne=number.edge(5);
x=GaussLobattoR(ne+2,0,1);
Re5=x(2:ne+1)-1;
Se5=Re5;
Te5=x(2:ne+1);
%% edge6
ne=number.edge(6);
x=GaussLobattoR(ne+2,0,1);
Re6=x(ne+1:-1:2);
Se6=x(2:ne+1)-1;
Te6=x(2:ne+1);
%% edge7
ne=number.edge(7);
x=GaussLobattoR(ne+2,0,1);
Re7=x(ne+1:-1:2);
Se7=Re7;
Te7=x(2:ne+1);
%% edge6
ne=number.edge(6);
x=GaussLobattoR(ne+2,0,1);
Re8=x(2:ne+1)-1;
Se8=x(ne+1:-1:2);
Te8=x(2:ne+1);
%% Vertex1-5
Rv=[-1;1;1;-1;0];
Sv=[-1;-1;1;1;0];
Tv=[0;0;0;0;1];
%%
R=[Rv;Re1;Re2;Re3;Re4;Re5;Re6;Re7;Re8;Rfq;Rf1;Rf2;Rf3;Rf4];
S=[Sv;Se1;Se2;Se3;Se4;Se5;Se6;Se7;Se8;Sfq;Sf1;Sf2;Sf3;Sf4];
T=[Tv;Te1;Te2;Te3;Te4;Te5;Te6;Te7;Te8;Tfq;Tf1;Tf2;Tf3;Tf4];
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

%% demo
% number.edge=[3,3,3,3,3,3,3,3];
% number.face={[10,10],3,3,3,3};
% number.H=[2,2,2];
% 
% [R,S,T]=PyramidFekete(number);
% X1=[-1,1,0];Y1=[-1,-1,0];Z1=[0,0,1];
% X2=[1,1,0];Y2=[-1,1,0];Z2=[0,0,1];
% X3=[1,-1,0];Y3=[1,1,0];Z3=[0,0,1];
% X4=[-1,-1,0];Y4=[1,-1,0];Z4=[0,0,1];
% X5=[-1,1,1,-1];Y5=[-1,-1,1,1];Z5=[0,0,0,0];
% X={X1,X2,X3,X4,X5};
% Y={Y1,Y2,Y3,Y4,Y5};
% Z={Z1,Z2,Z3,Z4,Z5};
% for i=1:5
% patch(X{i},Y{i},Z{i},[0.1*i,0.2*i,0.1*i]);
% alpha(0.5);
% hold on
% end
% scatter3(R,S,T,'filled','r');
% for i=1:length(S)
%    text(R(i)+0.02,S(i)+0.02, T(i)+0.02, num2str(i));
% end
% axis off
