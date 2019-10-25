function [G,Gr,Gs,Gt]=BaseTetra(x, y, z,number)
%% 求取四面体基的函数值和三个偏导数值
% number为结构体：储存四面体单元边、面、体的阶次
%  可缺省调用 [G]=BaseTetra(x,y,z,number)

ne=number.edge;n=number.H;
nf=[number.face{1},number.face{2},number.face{3},number.face{4}];

x=x(:); y=y(:); z=z(:);
if nargout==1
    C=[1-x-y-z, x, y, z];
    E1=TetraBasisE(x, y, z, ne(1), 2);
    E2=TetraBasisE(x, y, z, ne(2), 5);
    E3=TetraBasisE(x, y, z, ne(3), 4);
    E4=TetraBasisE(x, y, z, ne(4), 1);
    E5=TetraBasisE(x, y, z, ne(5), 3);
    E6=TetraBasisE(x, y, z, ne(6), 6);
    F1=TetraBasisF(x, y, z, nf(1), 1);
    F2=TetraBasisF(x, y, z, nf(2), 2);
    F3=TetraBasisF(x, y, z, nf(3), 4);
    F4=TetraBasisF(x, y, z, nf(4), 3);
    B=TetraBasisB(x, y, z, n, [1,1,1,1]);
    G=[C, E1, E2, E3, E4, E5, E6, F1, F2, F3, F4, B]; 
elseif nargout>1
    C=[1-x-y-z, x, y, z]; Cx=0*C; Cy=Cx; Cz=Cx;
    Cx(:,1)=-1; Cx(:,2)=1; 
    Cy(:,1)=-1; Cy(:,3)=1; 
    Cz(:,1)=-1; Cz(:,4)=1;     
    [E1, jace1]=TetraBasisE(x, y, z, ne(1), 2);
    [E2, jace2]=TetraBasisE(x, y, z, ne(2), 5);
    [E3, jace3]=TetraBasisE(x, y, z, ne(3), 4);
    [E4, jace4]=TetraBasisE(x, y, z, ne(4), 1);
    [E5, jace5]=TetraBasisE(x, y, z, ne(5), 3);
    [E6, jace6]=TetraBasisE(x, y, z, ne(6), 6);
    [F1, jacf1]=TetraBasisF(x, y, z, nf(1), 1);
    [F2, jacf2]=TetraBasisF(x, y, z, nf(2), 2);
    [F3, jacf3]=TetraBasisF(x, y, z, nf(3), 4);
    [F4, jacf4]=TetraBasisF(x, y, z, nf(4), 3);
    [B, jacb]=TetraBasisB(x, y, z, n, [1,1,1,1]);
    G=[C, E1, E2, E3, E4, E5, E6, F1, F2, F3, F4, B]; 
    Gx=[Cx, jace1{1}, jace2{1}, jace3{1}, jace4{1}, jace5{1}, jace6{1}, ...
            jacf1{1}, jacf2{1}, jacf3{1}, jacf4{1}, jacb{1}];
    Gy=[Cy, jace1{2}, jace2{2}, jace3{2}, jace4{2}, jace5{2}, jace6{2}, ...
            jacf1{2}, jacf2{2}, jacf3{2}, jacf4{2}, jacb{2}];
    Gz=[Cz, jace1{3}, jace2{3}, jace3{3}, jace4{3}, jace5{3}, jace6{3}, ...
            jacf1{3}, jacf2{3}, jacf3{3}, jacf4{3}, jacb{3}];  
    jac={Gx, Gy, Gz};
    Gr=jac{1};Gs=jac{2};Gt=jac{3};
end


%% demo - basis 1
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[2,2,2,2,2,2]+0; nf=[2,2,2,2]+0; n=3; N=15;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [G, jac]=TetraBasis(x, y, z, ne, nf, n);
% Gx=jac{1}; Gy=jac{2}; Gz=jac{3}; 
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% Gxi=Dx*G;
% Gyi=Dy*G;
% Gzi=Dz*G;
% 
% % Draw basis
% k=1;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, G(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Gy(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Gyi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');


%% demo - basis 2
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% ne=[2,2,2,2,2,2]+1; nf=[1,1,1,1]+1; n=2; 
% N=max([ne(:); nf(:)+3; n+4])+10;
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% [x, y, z]=meshgridtetrac(N);
% 
% % Get the basis on the unti tetrahedron
% [xi, yi, zi, trih]=tetranodes(ne, nf, n);
% Q=TetraBasis(xi, yi, zi, ne-2, nf, n);
% [G, jac]=TetraBasis(x, y, z, ne-2, nf, n);
% iQ=inv(Q); G=G*iQ; Gx=jac{1}*iQ; 
% Gy=jac{2}*iQ; Gz=jac{3}*iQ; 
% 
% % Draw basis
% k=1;
% Gk=reshape(G(:,k), N, N, N); 
% figure; hold on;
% tetrasurfg(x, y, z, Gk);
% plot3(xi, yi, zi, 'r.', 'MarkerSize', 20);
% axis equal;


%% Interpolation 1
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% ne=[2,2,2,2,2,2]+3; nf=[2,2,2,2]+3; n=4; N=35;
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne+2, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [xi, yi, zi, trih]=tetranodes(ne+2, nf, n);
% Q=TetraBasis(xi, yi, zi, ne, nf, n);
% [G, jac]=TetraBasis(x, y, z, ne, nf, n);
% iQ=inv(Q); G=G*iQ; Gx=jac{1}*iQ; 
% Gy=jac{2}*iQ; Gz=jac{3}*iQ; 
% 
% % Function values and derivatives
% fun=@(x,y,z) cos(x).*cos(y).*cos(z);
% funx=@(x,y,z) -sin(x).*cos(y).*cos(z);
% funy=@(x,y,z) -cos(x).*sin(y).*cos(z);
% funz=@(x,y,z) -cos(x).*cos(y).*sin(z);
% F0=fun(xi,yi,zi); 
% F=fun(x,y,z); Fx=funx(x,y,z); 
% Fy=funy(x,y,z); Fz=funz(x,y,z); 
% Fi=G*F0; Fxi=Gx*F0; Fyi=Gy*F0; Fzi=Gz*F0; 
% 
% % Draw basis
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Fz); 
% end
% plot3(xi, yi, zi, 'r.', 'MarkerSize', 20);
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Fzi-Fz); 
% end
% plot3(xi, yi, zi, 'r.', 'MarkerSize', 20);
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% title('DQH');


%% Interpolation 2
% % The number of edge nodes (ne), face nodes (nf), and
% %  the number of inside hierarchical basis (n)
% ne=[2,2,2,2,2,2]+7; nf=[2,2,2,2]+5; n=4; 
% N=max([ne(:); nf(:)+3; n+4])+10;
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get integration nodes and weights
% [x, y, z, tri]=GaussTetraH(ne, nf, n);
% 
% % Displacement functions on the tetrahedron
% Q=TetraBasis(x, y, z, ne-2, nf, n);
% [G, jac]=TetraBasis(x, y, z, ne-2, nf, n);
% iQ=inv(Q); G=G*iQ; Gx=jac{1}*iQ; 
% Gy=jac{2}*iQ; Gz=jac{3}*iQ; 
% 
% % Function values and derivatives
% fun=@(x,y,z) cos(x).*cos(y).*cos(z);
% funx=@(x,y,z) -sin(x).*cos(y).*cos(z);
% funy=@(x,y,z) -cos(x).*sin(y).*cos(z);
% funz=@(x,y,z) -cos(x).*cos(y).*sin(z);
% F0=fun(x, y, z); 
% F=fun(x, y, z); Fx=funx(x, y, z); 
% Fy=funy(x, y, z); Fz=funz(x, y, z); 
% Fi=G*F0; Fxi=Gx*F0; Fyi=Gy*F0; Fzi=Gz*F0; 
% 
% % Draw results
% figure; hold on;
% for k=1:4+n+2
%     trisurf(tri{k}, x, y, z, Fz); 
% end
% % plot3(x, y, z, 'r.', 'MarkerSize', 20);
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for k=1:4+n+2
%     trisurf(tri{k}, x, y, z, Fz-Fzi); 
% end
% % plot3(x, y, z, 'r.', 'MarkerSize', 20);
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% title('DQH');
end
%%
function [E, jac, ijk]=TetraBasisE(x, y, z, N, ed)

% TetraBasisE : Get edge functions on a unit tetrahedron and their derivatives
% 
% Calling Sequences:
%
%     E=TetraBasisE(x, y, z, N, ed)
%
%     [E, jac, ijk]=TetraBasisE(x, y, z, N, ed)
% 
% INPUTS:
% 
%     x, y, z   -    Vectors of parametric coordinates of a unit tetrahedron
%
%     n - The number of nodes on the unit tetrahedron
%
%     ed  -  Edge numbers from 1 to 6
%
% OUTPUT:
% 
%      E   -  Edge functions
%
%      jac - First derivatives of the edge functions
%
%      ijk - Index of i, j, k
%

if nargout==1
    switch ed
        case 1
            E=TetraBasisEk(x, y, z, N, [1,0,0,1]);
        case 2
            E=TetraBasisEk(x, y, z, N, [1,1,0,0]);
        case 3
            E=TetraBasisEi(x, y, z, N, [0,1,0,1]);
        case 4
            E=TetraBasisEk(x, y, z, N, [1,0,1,0]);
        case 5
            E=TetraBasisEi(x, y, z, N, [0,1,1,0]);
        case 6
            E=TetraBasisEj(x, y, z, N, [0,0,1,1]);
    end
elseif nargout>=2
    switch ed
        case 1
            [E, jac, ijk]=TetraBasisEk(x, y, z, N, [1,0,0,1]);
        case 2
            [E, jac, ijk]=TetraBasisEk(x, y, z, N, [1,1,0,0]);
        case 3
            [E, jac, ijk]=TetraBasisEi(x, y, z, N, [0,1,0,1]);
        case 4
            [E, jac, ijk]=TetraBasisEk(x, y, z, N, [1,0,1,0]);
        case 5
            [E, jac, ijk]=TetraBasisEi(x, y, z, N, [0,1,1,0]);
        case 6
            [E, jac, ijk]=TetraBasisEj(x, y, z, N, [0,0,1,1]);
    end
end
end
%% 
function [E, jac, ijk]=TetraBasisEi(x, y, z, N, pe)

% TetraBasisEi : Get tetrahedral basis functions and their derivatives on edges 3, 5
%
% INPUTS: 
%
%     N - The number of nodes on the edge
%
%     pe = [0,1,0,1] for edge 3
%
%     pe = [0,1,1,0] for edge 5
%

l=pe(2); m=pe(3); n=pe(4);
a=2*m+2*n+1; b=2*l; 

% Get the basis on the unti tetrahedron
E=zeros(length(x), N);
ijk=zeros(N, 3);
if nargout==1
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
    [Pi]=JacobiTetrai(N-1, a, b, x, y, z); 
elseif nargout>=2
    Ex=E; Ey=E; Ez=E; 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
    [Pi, jaci]=JacobiTetrai(N-1, a, b, x, y, z); 
end
t=1; j=0; k=0;
for i=0:N-1
    E(:,t)=Pp(:).*Pi(:,i+1);
    if nargout>=2
        Ex(:,t)=jacp{1}(:).*Pi(:,i+1)+Pp(:).*jaci{1}(:,i+1); 
        Ey(:,t)=jacp{2}(:).*Pi(:,i+1)+Pp(:).*jaci{2}(:,i+1); 
        Ez(:,t)=jacp{3}(:).*Pi(:,i+1)+Pp(:).*jaci{3}(:,i+1); 
    end
    ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
    t=t+1;
end
if nargout>=2
    jac={Ex, Ey, Ez};
end
end
%% 
function [E, jac, ijk]=TetraBasisEj(x, y, z, N, pe)

% TetraBasisEj : Get tetrahedral basis functions and their derivatives on edges 6
%
% INPUTS: 
%
%     N - The number of nodes on the edge
%
%     pe = [0,0,1,1] for edge 6
%

m=pe(3); n=pe(4);
c=2*n; d=2*m;

% Get the basis on the unti tetrahedron
E=zeros(length(x), N);
ijk=zeros(N, 3);
if nargout==1
    Pj=JacobiTetraj(N-1, c, d, x, y, z);
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
elseif nargout>=2
    Ex=E; Ey=E; Ez=E; 
    [Pj, jacj]=JacobiTetraj(N-1, c, d, x, y, z); 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
end
t=1; i=0; k=0;
for j=0:N-1
    E(:,t)=Pp(:).*Pj(:,j+1);
    if nargout>=2
        Ex(:,t)=jacp{1}(:).*Pj(:,j+1)+Pp(:).*jacj{1}(:,j+1); 
        Ey(:,t)=jacp{2}(:).*Pj(:,j+1)+Pp(:).*jacj{2}(:,j+1); 
        Ez(:,t)=jacp{3}(:).*Pj(:,j+1)+Pp(:).*jacj{3}(:,j+1); 
    end
    ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
    t=t+1;
end
if nargout>=2
    jac={Ex, Ey, Ez};
end
end
%%
function [E, jac, ijk]=TetraBasisEk(x, y, z, N, pe)

% TetraBasisEk : Get tetrahedral basis functions and their derivatives on edges 1, 2, 4
%
% INPUTS: 
%
%     N - The number of nodes on the edge
%
%     pe = [1,0,0,1] for edge 1
%
%     pe = [1,1,0,0] for edge 2
%
%     pe = [1,0,1,0] for edge 4
%

h=pe(1); l=pe(2); m=pe(3); n=pe(4);
e=2*h; f=2*l+2*m+2*n+2;

% Get the basis on the unti tetrahedron
E=zeros(length(x), N);
ijk=zeros(N, 3);
if nargout==1
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
    dJ=JacobiDerCellA(N-1, e, f, 2*(x+y+z)-1, 0);
    Pk=cell2mat(dJ(:,1)'); 
elseif nargout>=2
    Ex=E; Ey=E; Ez=E; 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
    dJ=JacobiDerCellA(N-1, e, f, 2*(x+y+z)-1, 1);
    Pk=cell2mat(dJ(:,1)'); 
    dPk=cell2mat(dJ(:,2)'); 
    jack={2*dPk, 2*dPk, 2*dPk};
end
t=1; i=0; j=0;
for k=0:N-1
    E(:,t)=Pp(:).*Pk(:,k+1);
    if nargout>=2
        Ex(:,t)=jacp{1}(:).*Pk(:,k+1)+Pp(:).*jack{1}(:,k+1);
        Ey(:,t)=jacp{2}(:).*Pk(:,k+1)+Pp(:).*jack{2}(:,k+1);
        Ez(:,t)=jacp{3}(:).*Pk(:,k+1)+Pp(:).*jack{3}(:,k+1);
    end
    ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
    t=t+1;
end
if nargout>=2
    jac={Ex, Ey, Ez};
end

%% demo
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[4,3,3,5,3,4]+2; nf=[2,2,2,2]+2; n=3; N=20;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [B, jac, ijk]=TetraBasisE(x, y, z, n, 6);
% Bx=jac{1}; By=jac{2}; Bz=jac{3}; 
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% Bxi=Dx*B;
% Byi=Dy*B;
% Bzi=Dz*B;
% 
% % Draw basis
% k=2;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, B(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bz(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bzi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');

end
%%
function [F, jac, ijk]=TetraBasisF(x, y, z, N, f)

% TetraBasisF : Get face functions on a unit tetrahedron and their derivatives
% 
% Calling Sequences:
%
%     F=TetraBasisF(x, y, z, N, f)
%
%     [F, jac, ijk]=TetraBasisF(x, y, z, N, f)
% 
% INPUTS:
% 
%     x, y, z   -    Vectors of parametric coordinates of a unit tetrahedron
%
%     n - The number of nodes on the unit tetrahedron
%
%     f  -  Face numbers from 1 to 4
%
% OUTPUT:
% 
%      F   -  Face functions
%
%      jac - First derivatives of the face functions
%
%      ijk - Index of i, j, k
%

if nargout==1
    switch f
        case 1
            F=TetraBasisFik(x, y, z, N, [1,1,1,0]);
        case 2
            F=TetraBasisFik(x, y, z, N, [1,1,0,1]);
        case 3
            F=TetraBasisFjk(x, y, z, N, [1,0,1,1]);
        case 4
            F=TetraBasisFij(x, y, z, N, [0,1,1,1]);
    end
elseif nargout>=2
    switch f
        case 1
            [F, jac, ijk]=TetraBasisFik(x, y, z, N, [1,1,1,0]);
        case 2
            [F, jac, ijk]=TetraBasisFik(x, y, z, N, [1,1,0,1]);
        case 3
            [F, jac, ijk]=TetraBasisFjk(x, y, z, N, [1,0,1,1]);
        case 4
            [F, jac, ijk]=TetraBasisFij(x, y, z, N, [0,1,1,1]);
    end
end
end


function [F, jac, ijk]=TetraBasisFij(x, y, z, N, pe)

% TetraBasisFij : Get tetrahedral basis functions and their derivatives on face 4
%
% INPUTS: 
%
%     N - The number of nodes on the face
%
%     pe = [0,1,1,1] for face 4
%

l=pe(2); m=pe(3); n=pe(4);
a=2*m+2*n+1; b=2*l; c=2*n; d=2*m;

% Get the basis on the unti tetrahedron
F=zeros(length(x), N*(N+1)/2);
ijk=zeros(N*(N+1)/2, 3);
if nargout==1
    Pj=JacobiTetraj(N-1, c, d, x, y, z);
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
elseif nargout>=2
    Fx=F; Fy=F; Fz=F; 
    [Pj, jacj]=JacobiTetraj(N-1, c, d, x, y, z); 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
end
t=1; k=0;
for i=0:N-1
    for j=0:N-1-i
        [Pi, jaci]=JacobiTetrai(i, a+2*j, b, x, y, z);
        F(:,t)=Pp(:).*Pi(:,end).*Pj(:,j+1);
        if nargout>=2
            Fx(:,t)=jacp{1}(:).*Pi(:,end).*Pj(:,j+1)+Pp(:).*jaci{1}(:,end).*Pj(:,j+1)+Pp(:).*Pi(:,end).*jacj{1}(:,j+1); 
            Fy(:,t)=jacp{2}(:).*Pi(:,end).*Pj(:,j+1)+Pp(:).*jaci{2}(:,end).*Pj(:,j+1)+Pp(:).*Pi(:,end).*jacj{2}(:,j+1); 
            Fz(:,t)=jacp{3}(:).*Pi(:,end).*Pj(:,j+1)+Pp(:).*jaci{3}(:,end).*Pj(:,j+1)+Pp(:).*Pi(:,end).*jacj{3}(:,j+1); 
        end
        ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
        t=t+1;
    end
end
if nargout>=2
    jac={Fx, Fy, Fz};
end
end

function [F, jac, ijk]=TetraBasisFik(x, y, z, N, pe)

% TetraBasisFik : Get tetrahedral basis functions and their derivatives on face 1, 2
%
% INPUTS: 
%
%     N - The number of nodes on the face
%
%     pe = [1,1,1,0] for face 1
%
%     pe = [1,1,0,1] for face 2
%

h=pe(1); l=pe(2); m=pe(3); n=pe(4);
a=2*m+2*n+1; b=2*l; 
e=2*h; f=2*l+2*m+2*n+2;

% Get the basis on the unti tetrahedron
F=zeros(length(x), N*(N+1)/2);
ijk=zeros(N*(N+1)/2, 3);
if nargout==1
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
elseif nargout>=2
    Fx=F; Fy=F; Fz=F; 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
end
t=1; j=0;
for i=0:N-1
    for k=0:N-1-i
        [Pi, jaci]=JacobiTetrai(i, a+2*j, b, x, y, z);
        [Pk, dPk]=JacobiRecDer(k, e, f+2*i+2*j, 2*(x+y+z)-1);
        jack={2*dPk, 2*dPk, 2*dPk};
        F(:,t)=Pp(:).*Pi(:,end).*Pk;
        if nargout>=2
            Fx(:,t)=jacp{1}(:).*Pi(:,end).*Pk+Pp(:).*jaci{1}(:,end).*Pk+Pp(:).*Pi(:,end).*jack{1};
            Fy(:,t)=jacp{2}(:).*Pi(:,end).*Pk+Pp(:).*jaci{2}(:,end).*Pk+Pp(:).*Pi(:,end).*jack{2};
            Fz(:,t)=jacp{3}(:).*Pi(:,end).*Pk+Pp(:).*jaci{3}(:,end).*Pk+Pp(:).*Pi(:,end).*jack{3};
        end
        ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
        t=t+1;
    end
end
if nargout>=2
    jac={Fx, Fy, Fz};
end
end

function [F, jac, ijk]=TetraBasisFjk(x, y, z, N, pe)

% TetraBasisFjk : Get tetrahedral basis functions and their derivatives on face 3
%
% INPUTS: 
%
%     N - The number of nodes on the face
%
%     pe = [1,0,1,1] for face 3
%

h=pe(1); l=pe(2); m=pe(3); n=pe(4);
c=2*n; d=2*m; e=2*h; f=2*l+2*m+2*n+2;

% Get the basis on the unti tetrahedron
F=zeros(length(x), N*(N+1)/2);
ijk=zeros(N*(N+1)/2, 3);
if nargout==1
    Pj=JacobiTetraj(N-1, c, d, x, y, z);
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
elseif nargout>=2
    Fx=F; Fy=F; Fz=F; 
    [Pj, jacj]=JacobiTetraj(N-1, c, d, x, y, z); 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
end
t=1; i=0;
for j=0:N-1
    for k=0:N-1-j
        [Pk, dPk]=JacobiRecDer(k, e, f+2*i+2*j, 2*(x+y+z)-1);
        jack={2*dPk, 2*dPk, 2*dPk};
        F(:,t)=Pp(:).*Pj(:,j+1).*Pk;
        if nargout>=2
            Fx(:,t)=jacp{1}(:).*Pj(:,j+1).*Pk+Pp(:).*jacj{1}(:,j+1).*Pk+Pp(:).*Pj(:,j+1).*jack{1}; 
            Fy(:,t)=jacp{2}(:).*Pj(:,j+1).*Pk+Pp(:).*jacj{2}(:,j+1).*Pk+Pp(:).*Pj(:,j+1).*jack{2}; 
            Fz(:,t)=jacp{3}(:).*Pj(:,j+1).*Pk+Pp(:).*jacj{3}(:,j+1).*Pk+Pp(:).*Pj(:,j+1).*jack{3}; 
        end
        ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
        t=t+1;
    end
end
if nargout>=2
    jac={Fx, Fy, Fz};
end

%% demo
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[4,3,3,5,3,4]+2; nf=[2,2,2,2]+2; n=3; N=20;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [B, jac, ijk]=TetraBasisF(x, y, z, n, 4);
% Bx=jac{1}; By=jac{2}; Bz=jac{3}; 
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% Bxi=Dx*B;
% Byi=Dy*B;
% Bzi=Dz*B;
% 
% % Draw basis
% k=5;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, B(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bz(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bzi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [B, jac, ijk]=TetraBasisB(x, y, z, N, pe)

% TetraBasisB : Get body functions on a unit tetrahedron and their derivatives
% 
% Calling Sequences:
%
%     B=TetraBasisN(x, y, z, N, pe)
%
%     [B, jac]=TetraBasisN(x, y, z, N, pe)
% 
% INPUTS:
% 
%     x, y, z   -    Vectors of parametric coordinates of a unit tetrahedron
%
%     n - The number of nodes on an edge of  the unit tetrahedron
%
%     pe =[1,1,1,1] -  Coefficients of the power of the polynomial of the unit tetrahedron
%
% OUTPUT:
% 
%      F   -  Face functions
%
%      jac - First derivatives of the Face functions
%

h=pe(1); l=pe(2); m=pe(3); n=pe(4);
a=2*m+2*n+1; b=2*l; c=2*n; d=2*m;
e=2*h; f=2*l+2*m+2*n+2;

% Get the basis on the unti tetrahedron
B=zeros(length(x), N*(N+1)*(N+2)/6);
ijk=zeros(N*(N+1)*(N+2)/6, 3);
if nargout==1
    Pj=JacobiTetraj(N-1, c, d, x, y, z);
    Pp=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4));
elseif nargout>=2
    Bx=B; By=B; Bz=B; 
    [Pj, jacj]=JacobiTetraj(N-1, c, d, x, y, z); 
    [Pp, jacp]=polytetrader(x, y, z, pe(1), pe(2), pe(3), pe(4)); 
end
t=1;
for i=0:N-1
    for j=0:N-1-i
        for k=0:N-1-i-j
            [Pi, jaci]=JacobiTetrai(i, a+2*j, b, x, y, z);
            [Pk, dPk]=JacobiRecDer(k, e, f+2*i+2*j, 2*(x+y+z)-1);
            jack={2*dPk, 2*dPk, 2*dPk};
            B(:,t)=Pp(:).*Pi(:,end).*Pj(:,j+1).*Pk;
            if nargout>=2
                Bx(:,t)=jacp{1}(:).*Pi(:,end).*Pj(:,j+1).*Pk+Pp(:).*jaci{1}(:,end).*Pj(:,j+1).*Pk+ ...
                    Pp(:).*Pi(:,end).*jacj{1}(:,j+1).*Pk+Pp(:).*Pi(:,end).*Pj(:,j+1).*jack{1};
                By(:,t)=jacp{2}(:).*Pi(:,end).*Pj(:,j+1).*Pk+Pp(:).*jaci{2}(:,end).*Pj(:,j+1).*Pk+ ...
                    Pp(:).*Pi(:,end).*jacj{2}(:,j+1).*Pk+Pp(:).*Pi(:,end).*Pj(:,j+1).*jack{2};
                Bz(:,t)=jacp{3}(:).*Pi(:,end).*Pj(:,j+1).*Pk+Pp(:).*jaci{3}(:,end).*Pj(:,j+1).*Pk+ ...
                    Pp(:).*Pi(:,end).*jacj{3}(:,j+1).*Pk+Pp(:).*Pi(:,end).*Pj(:,j+1).*jack{3};
            end
            ijk(t,1)=i; ijk(t,2)=j; ijk(t,3)=k; 
            t=t+1;
        end
    end
end
if nargout>=2
    jac={Bx, By, Bz};
end


%% demo
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[4,3,3,5,3,4]+2; nf=[2,2,2,2]+2; n=3; N=20;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [B, jac, ijk]=TetraBasisB(x, y, z, n, [1,1,1,1]);
% Bx=jac{1}; By=jac{2}; Bz=jac{3}; 
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% Bxi=Dx*B;
% Byi=Dy*B;
% Bzi=Dz*B;
% 
% % Draw basis
% k=2;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, B(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.3);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bz(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Bzi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [p, jac]=polytetrader(x, y, z, a, b, c, d)

% Get derivatives of polynomials of a unit tetrahedron: (1-x-y-z)^a*x^b*y^c*z^d
% 
% Calling Sequences:
%
%    [p, jac]=polytetrader(x, y, z, a, b, c, d)
% 
% INPUTS:
% 
%     x, y, z   -   Coordinates of parametric points.
%
%     a, b, c, d  - Power of the variables.
%
% OUTPUT:
% 
%     p   -   The value of the expression.
%
%     jac  - Evaluated first derivatives (Jacobian).
% 

% Get the value of (1-x-y-z)^a*x^b*y^c*z^d
if a<0 || b<0 || c<0 || d<0
    error('The order of the polynomial cannot be negative!');
end
p=polyeval([1-x-y-z, x, y, z], 1, [a,b,c,d]);

% First derivative
if nargout>=2
    px=polyeval([1-x-y-z, x, y, z], -a, [a-1,b,c,d])+polyeval([1-x-y-z, x, y, z], b, [a,b-1,c,d]);
    py=polyeval([1-x-y-z, x, y, z], -a, [a-1,b,c,d])+polyeval([1-x-y-z, x, y, z], c, [a,b,c-1,d]);
    pz=polyeval([1-x-y-z, x, y, z], -a, [a-1,b,c,d])+polyeval([1-x-y-z, x, y, z], d, [a,b,c,d-1]);
    jac={px, py, pz};
end

%% demo
% % Node numbers (n) and the order of the polynomials (a, b, c, d)
% N=17; a=2; b=2; c=2; d=2;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the value and derivatives of (1-x-y-z)^a*x^b*y^c*z^d
% [p, jac]=polytetrader(x, y, z, a, b, c, d);
% px=jac{1}; py=jac{2}; pz=jac{3}; 
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% pxi=Dx*p;
% pyi=Dy*p;
% pzi=Dz*p;
% 
% % Draw basis
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, p); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, py); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, pyi); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [Pi, jac]=JacobiTetrai(N, a, b, x, y, z)

% Get the Jacobi polynomials and their derivatives on a tetrahedron
% 
%  Calling sequence: 
% 
%    Pi=JacobiTetrai(N, a, b, x, y, z)
% 
%    [Pi, jac]=JacobiTetrai(N, a, b, x, y, z)
%   
%  Input :
%
%    N = The index of basis
%
%    a, b - Parameters of Jacobi polynomials
%
%    x, y, z - A vector of evaluation points on a unit tetrahedron
% 
%  Output :
% 
%    Pi - The Jacobi polynomials on a unit tetrahedron
% 
%    jac - The first derivatives of Jacobi polynomials on a unit tetrahedron
% 

% Coodinate transformation
s=x-y-z; t=x+y+z;

% Prepare matrices
n=length(s); 
Pi=zeros(n, N+1); 
if nargout>=2
    jac{1}=Pi; jac{2}=Pi; jac{3}=Pi; 
end

% Evaluate
if N==0
    Pi(:,1)=ones(n, 1);
elseif N==1
    P=t*(a-b)/2+s*(a+b+2)/2; 
    Pi(:,1)=ones(n, 1);
    Pi(:,2)=P;
    if nargout>=2
        dPs=(a+b+2)/2; 
        dPt=(a-b)/2; 
        jac{1}(:,2)=(dPt+dPs)*Pi(:,1);
        jac{2}(:,2)=(dPt-dPs)*Pi(:,1);
        jac{3}(:,2)=(dPt-dPs)*Pi(:,1);
    end
else
    P1=1; P2=t*(a-b)/2+s*(a+b+2)/2;
    Pi(:,1)=ones(n, 1); Pi(:,2)=P2;
    if nargout>=2
        dPs1=0; dPs2=(a+b+2)/2; 
        dPt1=0; dPt2=(a-b)/2; 
        jac{1}(:,2)=(dPt2+dPs2)*Pi(:,1);
        jac{2}(:,2)=(dPt2-dPs2)*Pi(:,1);
        jac{3}(:,2)=(dPt2-dPs2)*Pi(:,1);
    end
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b); 
        a2=(2*n+a+b+1)*(a^2-b^2); 
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2); 
        a4=2*(n+a)*(n+b)*(2*n+a+b+2); 
        P=((a2*t+a3*s).*P2-a4*t.^2.*P1)/a1; 
        Pi(:,n+2)=P; 
        if nargout>=2
            dPs=((a2*t+a3*s).*dPs2+a3*P2-a4*t.^2.*dPs1)/a1; 
            dPt=((a2*t+a3*s).*dPt2+a2*P2-a4*t.^2.*dPt1-2*a4*t.*P1)/a1; 
            jac{1}(:,n+2)=dPt+dPs; 
            jac{2}(:,n+2)=dPt-dPs; 
            jac{3}(:,n+2)=dPt-dPs;             
        end
        
        P1=P2; P2=P; 
        if nargout>=2
            dPs1=dPs2; dPs2=dPs; 
            dPt1=dPt2; dPt2=dPt; 
        end
    end
end


%% demo
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[4,3,3,5,3,4]+2; nf=[2,2,2,2]+2; n=4; N=15;
% a=3; b=2; i=3; j=4;
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [Pi, jac]=JacobiTetrai(i, a+2*j, b, x, y, z);
% dPix=jac{1}; dPiy=jac{2}; dPiz=jac{3};
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% dPixi=Dx*Pi;
% dPiyi=Dy*Pi;
% dPizi=Dz*Pi;
% 
% % Draw basis
% k=3;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Pi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, dPiy(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, dPiyi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [Pj, jac]=JacobiTetraj(N, a, b, x, y, z)

% Get the Jacobi polynomials and their derivatives on a tetrahedron
% 
%  Calling sequence: 
% 
%    Pj=JacobiTetraj(N, a, b, x, y, z)
% 
%    [Pj, jac]=JacobiTetraj(N, a, b, x, y, z)
%   
%  Input :
%
%    N = The index of basis
%
%    a, b - Parameters of Jacobi polynomials
%
%    x, y, z - A vector of evaluation points on a unit tetrahedron
% 
%  Output :
% 
%    Pj - The Jacobi polynomials on a unit tetrahedron
% 
%    jac - The first derivatives of Jacobi polynomials on a unit tetrahedron
% 

% Coodinate transformation
s=y-z; t=y+z;

% Prepare matrices
n=length(s); 
Pj=zeros(n, N+1); 
if nargout>=2
    jac{1}=Pj; jac{2}=Pj; jac{3}=Pj; 
end

% Evaluate
if N==0
    Pj(:,1)=ones(n, 1);
elseif N==1
    P=t*(a-b)/2+s*(a+b+2)/2; 
    Pj(:,1)=ones(n, 1);
    Pj(:,2)=P;
    if nargout>=2
        dPs=(a+b+2)/2; 
        dPt=(a-b)/2; 
        jac{1}(:,2)=0;
        jac{2}(:,2)=(dPt+dPs)*Pj(:,1);
        jac{3}(:,2)=(dPt-dPs)*Pj(:,1);
    end
else
    P1=1; P2=t*(a-b)/2+s*(a+b+2)/2;
    Pj(:,1)=ones(n, 1); Pj(:,2)=P2;
    if nargout>=2
        dPs1=0; dPs2=(a+b+2)/2; 
        dPt1=0; dPt2=(a-b)/2; 
        jac{1}(:,2)=0;
        jac{2}(:,2)=(dPt2+dPs2)*Pj(:,1);
        jac{3}(:,2)=(dPt2-dPs2)*Pj(:,1);
    end
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b); 
        a2=(2*n+a+b+1)*(a^2-b^2); 
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2); 
        a4=2*(n+a)*(n+b)*(2*n+a+b+2); 
        P=((a2*t+a3*s).*P2-a4*t.^2.*P1)/a1; 
        Pj(:,n+2)=P; 
        if nargout>=2
            dPs=((a2*t+a3*s).*dPs2+a3*P2-a4*t.^2.*dPs1)/a1; 
            dPt=((a2*t+a3*s).*dPt2+a2*P2-a4*t.^2.*dPt1-2*a4*t.*P1)/a1; 
            jac{1}(:,n+2)=0; 
            jac{2}(:,n+2)=dPt+dPs; 
            jac{3}(:,n+2)=dPt-dPs;             
        end
        
        P1=P2; P2=P; 
        if nargout>=2
            dPs1=dPs2; dPs2=dPs; 
            dPt1=dPt2; dPt2=dPt; 
        end
    end
end

%% demo
% % The number of edge nodes (ne), face nodes (nf), the number of
% %  inside hierarchical basis (n) and the node number for plot (N)
% % Parameters of Jacobi polynomials (a, b)
% ne=[4,3,3,5,3,4]+2; nf=[2,2,2,2]+2; n=4; N=15;
% a=2; b=2; i=5; 
% 
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
% 
% % Totla number of nodes
% [tn, tnf]=tetranumh(ne, nf, n);
% 
% % Get nodes
% u=GaussLobattoR(N, 0, 1);
% [x, y, z]=TetraLobatto(u);
% tri=tetradelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)']; 
% 
% % Get the basis on the unti tetrahedron
% [Pi, jac]=JacobiTetraj(i, a, b, x, y, z);
% dPix=jac{1}; dPiy=jac{2}; dPiz=jac{3};
% 
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% dPixi=Dx*Pi;
% dPiyi=Dy*Pi;
% dPizi=Dz*Pi;
% 
% % Draw basis
% k=3;
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, Pi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, dPiy(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
% 
% figure; hold on;
% for i=1:N+2
%     trisurf(tri{i}, x, y, z, dPiyi(:,k)); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function dJ=JacobiDerCellA(N, a, b, x, der)

% Get ALL of the N Jacobi polynomials and their derivatives
% 
%  Calling sequence: 
% 
%    dJ=JacobiDerCellA(N, a, b, x, der)
%   
%  Input :
%
%    N = the number of basis
%    a, b - parameters of Jacobi polynomials
%    x - a vector of evaluation points
%    der - the order of derivatives
% 
%  Output :
% 
%    dJ - the Jacobi polynomials and their derivatives
% 

dJ=cell(N+1, der+1); 

if N==0
    dJ{1}=ones(size(x));
    for t=2:der+1
        dJ{t}=0*dJ{1};
    end
elseif N==1
    dJ(1,:)=JacobiDerCellA(0, a, b, x, der);
    dJ{2,1}=0.5*(a-b+(a+b+2)*x);
    if der>0
        dJ{2,2}=0.5*(a+b+2)*ones(size(x));
    end
    for t=3:der+1
        dJ{2,t}=0*dJ{2,1};
    end
elseif N>1
    pp=JacobiDerCellA(1, a, b, x, der);
    dJ([1,2],:)=pp(:,1:der+1);
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b); 
        a2=(2*n+a+b+1)*(a^2-b^2); 
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2); 
        a4=2*(n+a)*(n+b)*(2*n+a+b+2); 
        dJ{n+2,1}=((a2+a3*x).*dJ{n+1,1}-a4*dJ{n,1})/a1; 
        for t=2:der+1
            dJ{n+2,t}=((a2+a3*x).*dJ{n+1,t}+(t-1)*a3*dJ{n+1,t-1}-a4*dJ{n,t})/a1;   
        end
    end
end


%% demo
% % The number of points (n), the order of basis (m), 
% %    and the order of derivatives (der)
% n=100; m=3; dN=4; a=3; b=3; der=2;
% 
% % Get Gauss-Birkhoff points
% x=LobattoChebyshev(-1,1,n);  
% A=Weighting(x); 
% 
% dJa=JacobiDerCellA(m-dN, a, b, x, der);
% dJa=JacobiDerCellH(dJa, dN, a, b, x);
% dJ=dJa(m+1,:);
% y=dJ{1}; dy=dJ{2}; ddy=dJ{3}; 
% 
% % Symbolic Jacobi polynomials
% syms s
% Jh=jacobiP(m, a, b, s);    
% rt=roots(sym2poly(Jh));
% 
% % st=dJacobiRt(k+1, a-1, b-1);
% st=JacobiRt(m,a,b);
% figure; hold on;
% plot(x,y); 
% plot(rt,0*rt,'ko');
% plot(st,0*st,'r*');
% plot([-1, 1], [0, 0]);
% title('J');
% 
% st=JacobiRt(m-1, a+1, b+1);
% figure; hold on
% plot(x,dy); 
% plot(st,0*st,'ro');
% plot([-1, 1], [0, 0]);
% title('dJ');
% 
% st=JacobiRt(m-2, a+2, b+2);
% figure; hold on
% plot(x,ddy); 
% plot(st,0*st,'ro');
% plot([-1, 1], [0, 0]);
% title('ddJ');
end
function [y, dy, ddy, dddy]=JacobiRecDer(N, a, b, x)

% Get the Jacobi polynomials and their derivatives
% 
%  Calling sequence: 
% 
%    y=JacobiRecDer(N, a, b, x)
% 
%    [y, dy]=JacobiRecDer(N, a, b, x)
% 
%    [y, dy, ddy]=JacobiRecDer(N, a, b, x)
% 
%    [y, dy, ddy, dddy]=JacobiRecDer(N, a, b, x)
%   
%  Input :
%
%    N = the number of basis
%    a, b - parameters of Jacobi polynomials
%    x - a vector of nodes
% 
%  Output :
% 
%    y - the Jacobi polynomials
%    dy, ddy, dddy - the first, second and third 
%          order derivatives of Jacobi polynomials
% 

if N<0
    y=zeros(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
elseif N==0
    y=ones(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
elseif N==1
    y=0.5*(a-b+(a+b+2)*x); dy=0.5*(a+b+2); ddy=0*y; dddy=0*y;
else
    y1=1; y2=0.5*(a-b+(a+b+2)*x);
    dy1=0; dy2=0.5*(a+b+2);
    ddy1=0; ddy2=0;
    dddy1=0; dddy2=0;
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b);
        a2=(2*n+a+b+1)*(a^2-b^2);
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
        a4=2*(n+a)*(n+b)*(2*n+a+b+2);
        y=((a2+a3*x).*y2-a4*y1)/a1;
        if nargout>1
            dy=((a2+a3*x).*dy2+a3*y2-a4*dy1)/a1; 
        end
        if nargout>2
            ddy=((a2+a3*x).*ddy2+2*a3*dy2-a4*ddy1)/a1;          
        end
        if nargout>3
            dddy=((a2+a3*x).*dddy2+3*a3*ddy2-a4*dddy1)/a1;   
        end
        
        y1=y2; y2=y;
        if nargout>1
            dy1=dy2; dy2=dy;
        end
        if nargout>2
            ddy1=ddy2; ddy2=ddy;
        end
        if nargout>3
            dddy1=dddy2; dddy2=dddy;
        end
    end
end
end
function z=polyeval(v, c, p)

% Evaluate c*x^a*y^b..., x=v(1,:), y=v(2,:), ...; a=p(1), b=p(2), ...
% 
% Calling Sequences:
%
%    z=polyeval(v, c, p)
% 
% INPUTS:
% 
%     v   -   A matrix of dim*n. Each vector corresponds to a variable.
%
%     c  - The coefficient of the polynomial.
%
%     p -  The power for each variable.
%
% OUTPUT:
% 
%     z   -   The value of the polynomial expression.
% 
% Discriptions:
%
%     See also polytrider.
%

[m, k]=size(v); n=length(p);
if m==n
    if c==0 || ~isempty(find(p<0, 1))
        z=zeros(1,k);
    else
        z=c;
        for i=1:n
            z=z.*v(i,:).^p(i);
        end
    end
elseif k==n
    if c==0 || ~isempty(find(p<0, 1))
        z=zeros(m,1);
    else
        z=c;
        for i=1:n
            z=z.*v(:,i).^p(i);
        end
    end
else
    error('The input for variables are not correct!');
end
end

