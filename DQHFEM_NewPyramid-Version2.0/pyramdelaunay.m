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


%% demo
% % The number of nodes defined on an edge (n)
% n=6; 
% 
% % Generate non-equally spaces nodes on a unit right pyramid
% [R, S, T, id]=PyramLobatto(n);
% quad=pyramdelaunay(n);
% 
% % Transfor the pyramid to be equilateral
% R=2*R-1+T;
% S=2*S-1+T;
% T=sqrt(2)*T;
% 
% % Plot
% figure; hold on;
% for i=1:n+3
%     trimesh(quad{i}, R, S, T);
% end
% % plot3(R, S, T, 'r.', 'MarkerSize', 20);
% for i=1:length(R)
%     text(R(i), S(i), T(i), num2str(i));
% end
% view(3);
% axis equal;




