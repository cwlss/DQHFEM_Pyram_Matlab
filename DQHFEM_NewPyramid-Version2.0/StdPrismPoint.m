function [X,Y,Z,tri]=StdPrismPoint(n)
%% 创建标准金字塔体离散点。
% L=n*(n+1)*(2*n+1)/6.0;
% zz=GaussLobattoR(n, 0, 1);
x=0:1/(n-1):1;
y=0:1/(n-1):1;
z=0:1/(n-1):1;
N=n*(n+1)/2;
xx=zeros(N,1);yy=xx;
ij=1;
for j=1:n
    for i=1:n-j+1
        xx(ij)=x(i);
        yy(ij)=y(j);
        ij=ij+1;
    end
end
X=[];Y=[];Z=[];
for k=1:n
    zz=ones(size(xx))*z(k);
    X=[X;xx];
    Y=[Y;yy];
    Z=[Z;zz];
end
%% 创建tridurf索引
tri=PrismDelaunay(n);
end

function tri=PrismDelaunay(n)
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
tri=cell(n+3,1);
% The Trangles of each layer
for i=1:n
    trii=tridelaunay(n);
    sn=(i-1)*n*(n+1)/2;
    tri{i}=trii+sn;
end

% The face 1
tri1=zeros((n-1)^2,1);
L=n*(n+1)/2;
ij=1;
for j=1:n-1
    for i=1:n-1
        v1=(j-1)*L+i;
        v2=(j-1)*L+i+1;
        v3=j*L+i+1;
        v4=j*L+i;
        tri1(ij,1)=v1;tri1(ij,2)=v2;tri1(ij,3)=v3;
        ij=ij+1;
        tri1(ij,1)=v4;tri1(ij,2)=v3;tri1(ij,3)=v1;
        ij=ij+1;
    end
end
tri{n+1}=tri1;

% face 2
tri2=zeros((n-1)^2,1);
ij=1;
L=n*(n+1)/2;
for j=1:n-1
    for i=1:n-1
        v1=(j-1)*L+(2*n-i+1)*i/2;
        v2=(j-1)*L+(2*n-i)*(i+1)/2;
        v3=j*L+(2*n-i)*(i+1)/2;
        v4=j*L+(2*n-i+1)*i/2;
        tri2(ij,1)=v1;tri2(ij,2)=v2;tri2(ij,3)=v3;
        ij=ij+1;
        tri2(ij,1)=v4;tri2(ij,2)=v3;tri2(ij,3)=v1;
        ij=ij+1;
    end
end
tri{n+2}=tri2;

% The triangular face 3
tri3=zeros((n-1)^2,1);
ij=1;
L=n*(n+1)/2;
for j=1:n-1
    for i=1:n-1
        v1=(j-1)*L+(2*n-i+2)*(i-1)/2+1;
        v2=(j-1)*L+(2*n-i+1)*i/2+1;
        v3=j*L+(2*n-i+1)*i/2+1;
        v4=j*L+(2*n-i+2)*(i-1)/2+1;
        tri3(ij,1)=v1;tri3(ij,2)=v2;tri3(ij,3)=v3;
        ij=ij+1;
        tri3(ij,1)=v4;tri3(ij,2)=v3;tri3(ij,3)=v1;
        ij=ij+1;
    end
end
tri{n+3}=tri3;

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










