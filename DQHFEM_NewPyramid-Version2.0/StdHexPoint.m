function [X,Y,Z,tri]=StdHexPoint(n)
%% 创建标准金字塔体离散点（均匀分布离散点，仅用于离散位移场的可视化，不参与插值过程）。
x=0:1/(n-1):1;
y=0:1/(n-1):1;
z=0:1/(n-1):1;
X=[];Y=[];Z=[];
[y, x]=meshgrid(y, x);
x=x(:); y=y(:);
for i=1:n
    zz=ones(n^2,1)*z(i);
    X=[X;x];
    Y=[Y;y];
    Z=[Z;zz];
end


%% 创建tridurf索引
tri=HexDelaunay(n);
end

function tri=HexDelaunay(n)

% pyramdelaunay: Delaunay triangulation of gridded nodes on a unit pyramid
%
% Calling Sequences:
% 
%       quad=HexDelaunay(n)
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
tri=cell(n+4,1);
for i=1:n
    trii=quadelaunay(n, n);
    sn=n*n*(i-1);
    tri{i}=trii+sn;
end
% side face 1 
tri1=zeros(2*(n-1)^2,3);
ij=1;
for i=1:n-1
    for j=1:n-1
        r1=(j-1)*n^2+i;r2=(j-1)*n^2+i+1;r3=j*n^2+i+1;
        tri1(ij,1)=r1;tri1(ij,2)=r2;tri1(ij,3)=r3;
        ij=ij+1;
        s1=j*n^2+i;s2=j*n^2+i+1;s3=(j-1)*n^2+i;
        tri1(ij,1)=s1;tri1(ij,2)=s2;tri1(ij,3)=s3;
        ij=ij+1;
    end
end
tri{n+1}=tri1;
% side face 2 
tri2=zeros(2*(n-1)^2,3);
ij=1;
for i=1:n-1
    for j=1:n-1
        r1=(j-1)*n^2+i*n;r2=(j-1)*n^2+(i+1)*n;r3=j*n^2+(i+1)*n;
        tri2(ij,1)=r1;tri2(ij,2)=r2;tri2(ij,3)=r3;
        ij=ij+1;
        s1=j*n^2+i*n;s2=r3;s3=r1;
        tri2(ij,1)=s1;tri2(ij,2)=s2;tri2(ij,3)=s3;
        ij=ij+1;
    end
end
tri{n+2}=tri2;
% side face 3
tri3=zeros(2*(n-1)^2,3);
ij=1;
for i=1:n-1
    for j=1:n-1
        v1=(j-1)*n^2+n*(n-1)+i;
        v2=(j-1)*n^2+n*(n-1)+i+1;
        v3=j*n^2+n*(n-1)+i+1;
        v4=j*n^2+n*(n-1)+i;
        tri3(ij,1)=v1;tri3(ij,2)=v2;tri3(ij,3)=v3;
        ij=ij+1;
        tri3(ij,1)=v4;tri3(ij,2)=v3;tri3(ij,3)=v1;
        ij=ij+1;
    end
end
tri{n+3}=tri3;
% side face4
tri4=zeros(2*(n-1)^2,3);
ij=1;
for i=1:n-1
    for j=1:n-1
        v1=(j-1)*n^2+n*(i-1)+1;
        v2=(j-1)*n^2+n*i+1;
        v3=j*n^2+n*i+1;
        v4=j*n^2+n*(i-1)+1;
        tri4(ij,1)=v1;tri4(ij,2)=v2;tri4(ij,3)=v3;
        ij=ij+1;
        tri4(ij,1)=v4;tri4(ij,2)=v3;tri4(ij,3)=v1;
        ij=ij+1;
    end
end
tri{n+4}=tri4;
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

%% demo
% n=10;
% [X,Y,Z,tri]=StdHexPoint(n);
% scatter3(X,Y,Z);
% figure;
% for i=1:n+4
%     trisurf(tri{i},X,Y,Z);
%     hold on
% end









