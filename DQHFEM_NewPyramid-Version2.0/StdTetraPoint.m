function [X,Y,Z,tri]=StdTetraPoint(n)
%% 创建标准四面体离散点。
% L=n*(n+1)*(2*n+1)/6.0;
% zz=GaussLobattoR(n, 0, 1);
z=0:1/(n-1):1;
L=n*(n+1)*(n+2)/6;
X=zeros(L,1);Y=X;Z=X;
ij=1;
for k=1:n-1
    p=n-k+1;
    x=0:(1-z(k))/(p-1):1-z(k);
    y=0:(1-z(k))/(p-1):1-z(k);
    for j=1:p
        for i=1:p-j+1
            X(ij)=x(i);
            Y(ij)=y(j);
            Z(ij)=z(k);
            ij=ij+1;
        end
    end
end
X(L)=0;Y(L)=0;Z(L)=1;
%% 创建tridurf索引
tri=TetraDelaunay(n);
end


function tri=TetraDelaunay(n)


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

% The triangles of each layer
L=n*(n+1)*(n+2)/6;
tri=cell(n+2,1);
sn=0;
for i=1:n-1
    p=n-i+1;
    trii=tridelaunay(p);
    tri{i}=trii+sn;
    sn=sn+p*(p+1)/2;
end
% The triangular face 1
tri1=zeros((n-1)^2,3);
sn=0;snn=0;
ij=1;
for j=1:n-2
    snn=snn+(n-j+2)*(n-j+1)/2;
    for i=1:n-j-1
        v1=sn+i;
        v2=sn+i+1;
        v3=snn+i+1;
        v4=snn+i;
        tri1(ij,1)=v1;tri1(ij,2)=v2;tri1(ij,3)=v4;
        ij=ij+1;
        tri1(ij,1)=v4;tri1(ij,2)=v3;tri1(ij,3)=v2;
        ij=ij+1;
    end
    tri1(ij,1)=v1+1;tri1(ij,2)=v2+1;tri1(ij,3)=v4+1;
    ij=ij+1;
    sn=sn+(n-j+2)*(n-j+1)/2;
end
tri1((n-1)^2,1)=sn+1;
tri1((n-1)^2,2)=sn+2;
tri1((n-1)^2,3)=sn+4;
tri{n}=tri1;
% The triangular face 2
tri2=zeros((n-1)^2,3);
ij=1;
sn=0;
for j=1:n-2
    ln=0;lnn=0;
    for i=1:n-j-1
        v1=sn+ln+1;
        v2=sn+ln+1+n-j-i+2;
        v4=sn+lnn+(n-j+1)*(n-j+2)/2+1;
        v3=sn+lnn+(n-j+1)*(n-j+2)/2+n-j-i+1+1;
        tri2(ij,1)=v1;tri2(ij,2)=v2;tri2(ij,3)=v4;
        ij=ij+1;
        tri2(ij,1)=v4;tri2(ij,2)=v3;tri2(ij,3)=v2;
        ij=ij+1;
        ln=ln+n-j-i+2;
        lnn=lnn+n-j-i+1;
    end
    sn=sn+(n-j+1)*(n-j+2)/2;
    tri2(ij,1)=v1+3;tri2(ij,2)=v2+2;tri2(ij,3)=v3;
    ij=ij+1;
end
tri2((n-1)^2,1)=L-3;
tri2((n-1)^2,2)=L-1;
tri2((n-1)^2,3)=L;
tri{n+1}=tri2;
% The triangular face 3
tri3=zeros((n-1)^2,3);
ij=1;
sn=0;
for j=1:n-2
    ln=0;
    for i=1:n-j-1
        ln=ln+n-j+1-i+1;
        lnn=ln+n-j-i+1;
        v1=sn+ln;
        v2=sn+lnn;
        v4=sn+(n-j+1)*(n-j+2)/2+ln-i;
        v3=v4+n-j-i;
        tri3(ij,1)=v1;tri3(ij,2)=v2;tri3(ij,3)=v4;
        ij=ij+1;
        tri3(ij,1)=v4;tri3(ij,2)=v3;tri3(ij,3)=v2;
        ij=ij+1;
    end
    sn=sn+(n-j+1)*(n-j+2)/2;
    tri3(ij,1)=v1+2;tri3(ij,2)=v2+1;tri3(ij,3)=v3;
    ij=ij+1;
end
tri3((n-1)^2,1)=L-2;
tri3((n-1)^2,2)=L-1;
tri3((n-1)^2,3)=L;
tri{n+2}=tri3;
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










