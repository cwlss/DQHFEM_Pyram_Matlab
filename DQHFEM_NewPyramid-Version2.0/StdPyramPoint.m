function [X,Y,Z,tri]=StdPyramPoint(n)
%% 创建标准金字塔体离散点。
% L=n*(n+1)*(2*n+1)/6.0;
% zz=GaussLobattoR(n, 0, 1);
zz=0:1/(n-1):1;
X=[];Y=[];Z=[];
for i=1:n-2
    p=n-i+1;
    %x=GaussLobattoR(p, zz(i)-1, 1-zz(i));
    x=zz(i)-1: 2*(1-zz(i))/(p-1):1-zz(i);
    %y=GaussLobattoR(p, zz(i)-1, 1-zz(i));
    y=zz(i)-1: 2*(1-zz(i))/(p-1):1-zz(i);
    [y, x]=meshgrid(y, x);
    x=x(:); y=y(:);
    z=ones(size(x))*zz(i);
    X=[X;x];
    Y=[Y;y];
    Z=[Z;z];
end
r=[zz(n-1)-1;1-zz(n-1);zz(n-1)-1;1-zz(n-1);0];
s=[zz(n-1)-1;zz(n-1)-1;1-zz(n-1);1-zz(n-1);0];
t=[zz(n-1);zz(n-1);zz(n-1);zz(n-1);1];
X=[X;r];
Y=[Y;s];
Z=[Z;t];
%% 创建tridurf索引
tri=pyramdelaunay(n);
end

function tri=pyramdelaunay(n)

% pyramdelaunay: Delaunay triangulation of gridded nodes on a unit pyramid
%
% Calling Sequences:
% 
%       quad=pyramdelaunay(n)
% 
% INPUTS:
% 
%       n - The number of nodes on edges. All four edges have the same number.      
%
% OUTPUT:
% 
%     tri - Cell array of triangulation connectivity list. The first n-1
%            are layers of the tetrahedron. From n to n+2 are the 3 
%            surfaces of the tetrahedron without the bottum layer.
%
% Discritopn:
%  
%    See also: PyramLobatto
% 

% The quadrangles of each layer
for i=1:n-1
    p=n-i+1;
    trii=quadelaunay(p, p);
    sn=n*(n+1)*(2*n+1)/6.0-p*(p+1)*(2*p+1)/6.0;
    tri{i}=trii+sn;
end
tri0=tridelaunay(n);
tn=sn+5;

% The triangular face 1
num=zeros(n*(n+1)/2, 1);
t=1;
for i=1:n-1
    p=n-i+1;
    sn=n*(n+1)*(2*n+1)/6.0-p*(p+1)*(2*p+1)/6.0;
    for j=1:p
        num(t)=sn+j;
        t=t+1;
    end
end
num(t)=tn;
trii=tri0;
for i=1:size(tri0,1)
    for j=1:3
        trii(i,j)=num(tri0(i,j));
    end
end
tri{n}=trii;

% The triangular face 2
t=1;
for i=1:n-1
    p=n-i+1;
    sn=n*(n+1)*(2*n+1)/6.0-p*(p+1)*(2*p+1)/6.0;
    for j=1:p
        num(t)=sn+j*p;
        t=t+1;
    end
end
num(t)=tn;
for i=1:size(tri0,1)
    for j=1:3
        trii(i,j)=num(tri0(i,j));
    end
end
tri{n+1}=trii;

% The triangular face 3
t=1;
for i=1:n-1
    p=n-i+1;
    sn=n*(n+1)*(2*n+1)/6.0-p*(p+1)*(2*p+1)/6.0;
    for j=1:p
        num(t)=sn+(p-1)*p+j;
        t=t+1;
    end
end
num(t)=tn;
for i=1:size(tri0,1)
    for j=1:3
        trii(i,j)=num(tri0(i,j));
    end
end
tri{n+2}=trii;

% The triangular face 4
t=1;
for i=1:n-1
    p=n-i+1;
    sn=n*(n+1)*(2*n+1)/6.0-p*(p+1)*(2*p+1)/6.0;
    for j=1:p
        num(t)=sn+(j-1)*p+1;
        t=t+1;
    end
end
num(t)=tn;
for i=1:size(tri0,1)
    for j=1:3
        trii(i,j)=num(tri0(i,j));
    end
end
tri{n+3}=trii;
end
function tri=quadelaunay(m, n)

% quadelaunay: Delaunay triangulation of nodes on a unit quadrangle.
%
% Calling sequence: 
% 
%      tri=quadelaunay(m, n)
% 
%  Input :
%
%    m, n - The number of nodes inside the element 
%               on s- and t-directions
% 
%  Output :
% 
%    tri - Triangulation connectivity list, specified as an m-by-dim matrix, 
%            where m is the number of triangles, and dim is the number of 
%            vertices per triangle. Each element in tri is a Vertex ID. 
%            Each row of tri contains the vertex IDs that define a triangle.
%

if m==0 || n==0
    tri=[];
else
    t=1;
    for j=1:n-1
        q1=(j-1)*m+1:j*m;
        q2=j*m+1:(j+1)*m;
        for i=1:m-1
            tri(t,1)=q1(i); tri(t,2)=q1(i+1); tri(t,3)=q2(i+1);  t=t+1;
            tri(t,1)=q2(i); tri(t,2)=q2(i+1); tri(t,3)=q1(i);  t=t+1;
        end
    end
end
%% demo
% m=10; n=12;
% x=linspace(0, 1, m); 
% y=linspace(0, 1, n); 
% [y, x]=meshgrid(y, x);
% x=x(:); y=y(:);
% tri=quadelaunay(m, n);
% figure; hold on;
% triplot(tri, x, y);
end
function tri=tridelaunay(n)

% tridelaunay: Delaunay triangulation of gridded nodes on a unit triangle.
%
% Calling Sequences:
% 
%       tri=tridelaunay(n)
% 
% INPUTS:
% 
%       n - The number of nodes on edges. All three edges have the same number.      
%
% OUTPUT:
% 
%     tri - Triangulation connectivity list, specified as an m-by-dim matrix, 
%            where m is the number of triangles, and dim is the number of 
%            vertices per triangle. Each element in tri is a Vertex ID. 
%            Each row of tri contains the vertex IDs that define a triangle.
%
% Discritopn:
%  
%    See also: TrigNodeVect
% 

% Get the matrix
p=1; Num=nan(n); 
tri=zeros((n-1)^2,3);
for j=1:n
    for i=1:(n-j)+1
        Num(j,i)=p;
        p=p+1;
    end
end
p=1;
for j=1:n-1
    for i=1:(n-j)-1
        tri(p,1)=Num(j,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i); 
        p=p+1;
        tri(p,1)=Num(j+1,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i+1); 
        p=p+1;
    end
    if isempty(i)
        i=1;
    else
        i=i+1;
    end
    tri(p,1)=Num(j,i); tri(p,2)=Num(j,i+1); tri(p,3)=Num(j+1,i); 
    p=p+1;
end
end

%% demo
% n=6;
% [X,Y,Z]=StdPyramidPoint(n);
% scatter3(X,Y,Z);










