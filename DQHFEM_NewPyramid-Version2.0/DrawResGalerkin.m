
function DrawResGalerkin(n)
% n=18; 
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],n-2,n-2,n-2,n-2};
number.H=n-2; M=n+2;

% Integration nodes
na=4; nx=M+na; ny=M+na; nz=M+na;
[x, y, z, c]=meshgridpyramc(nx, ny, nz);

% Get the basis on the unti tetrahedron
 [G, Gx, Gy, Gz]=BasePyram(x(:), y(:), z(:), number);

% Function values and derivatives
%% Function Test1
fun=@(x,y,z) cos(x).*cos(y).*cos(z);
funx=@(x,y,z) -sin(x).*cos(y).*cos(z);
funy=@(x,y,z) -cos(x).*sin(y).*cos(z);
funz=@(x,y,z) -cos(x).*cos(y).*sin(z);
%% Function Test2
% fun=@(x,y,z) (x.^3).*(y.^3).*(z.^8)+(x.^2).*(y)+2*(x)+(y)+1;
% funx=@(x,y,z) (3*x.^2).*(y.^3).*(z.^8)+(2*x).*(y)+2;
% funy=@(x,y,z) (x.^3).*(3*y.^2).*(z.^8)+(x.^2)+1;
% funz=@(x,y,z) (x.^3).*(y.^3).*(8*z.^7);
%% Function Test3
% fun=@(x,y,z) (x.^3)+(z.^8);
% funx=@(x,y,z) (3*x.^2);
% funy=@(x,y,z) 0;
% funz=@(x,y,z) (8*z.^7);
%% Function Test4
% fun=@(x,y,z) x.*y./(1.1-z);
% funx=@(x,y,z) y./(1.1-z);
% funy=@(x,y,z) x./(1.1-z);
% funz=@(x,y,z) x.*y./((1.1-z).^2);
%% Function Test5
% fun=@(x,y,z) x+y+z.^7+100.*z;
% funx=@(x,y,z) 0*x+1;
% funy=@(x,y,z) 0*x+1;
% funz=@(x,y,z) 7*z.^6+100;

F=fun(x,y,z); Fx=funx(x,y,z); 
Fy=funy(x,y,z); Fz=funz(x,y,z); 

% Galerkin interpolation

CJ=diag(c(:)); 
Gk=(G'*CJ*G)\G'*CJ; 
w=Gk*F(:); 

% % Plot results
figure; hexasurfg(x, y, z, F, 2); 
title('F'); axis off
figure; hexasurfg(x, y, z, Fx, 2); 
title('Fx - Exact'); axis off
figure; hexasurfg(x, y, z, Fy, 2); 
title('Fy - Exact'); axis off
figure; hexasurfg(x, y, z, Fz, 2); 
title('Fz - Exact'); axis off

Fi=reshape(G*w, nx, ny, nz);
figure; hexasurfg(x, y, z, F-Fi, 2); 
title('Res-F-Galerkin');
axis off

Fxi=reshape(Gx*w, nx, ny, nz);
figure; hexasurfg(x, y, z, Fx-Fxi, 2); 
title('Res-Fx-Galerkin'); 
axis off

Fyi=reshape(Gy*w, nx, ny, nz);
figure; hexasurfg(x, y, z, Fy-Fyi, 2); 
title('Res-Fy-Galerkin');
axis off

Fzi=reshape(Gz*w, nx, ny, nz);
figure; hexasurfg(x, y, z, Fz-Fzi, 2); 
title('Res-Fz-Galerkin');
axis off
end

function [R, S, T, C]=meshgridpyramc(tt, cc, K)

% Replicates the grid vectors rgv, sgv, tgv and the associated 
%     integration weights to produce the coordinates of 
%     a pyramid grid (R, S, T).
% 
% Calling Sequences:
%
%     [R, S, T, C]=meshgridpyramidc(tt, cc)
%
% INPUTS:
% 
%      tt   -  If tt is a cell array, tt = {tu, tv, tw} of the parametric coordinates.
%            All tu, tv and tw should be defined in [0, 1].
%              If both inputs are nunvers, tt should be the number intergration 
%              nodes on an edge of (r, s) plane, and cc should be the number 
%              intergration nodes on t direction.
%
%     cc  - integration weights in one dimensional for u, v and w direction 
%             respectively, cc is a cell {cu, cv, cw}.
%
%     K - Int number of nodes on t direction.
%
% OUTPUT:
% 
%     R    :   a matrix of u coordinates for integration on a pyramid
%     S    :   a matrix of v coordinates for integration on a pyramid
%     T    :   a matrix of w coordinates for integration on a pyramid
%     C    :   a matrix of integration weights for pyramid
% 

if nargin==3
    if isnumeric(tt) && isnumeric(cc)
        M=tt; N=cc;
        [r, Cr]=GaussLobattoR(M, -1, 1);
        [s, Cs]=GaussLobattoR(N, -1, 1);
        [t, Ct]=GaussLobattoR(K, 0, 1);
        tt={r, s, t}; cc={Cr, Cs, Ct};
    end
end

[s, r, t]=meshgrid(tt{2}, tt{1}, tt{3}); 
[Cs, Cr, Ct]=meshgrid(cc{2}, cc{1}, cc{3}); 

% Get nodes on a unit pyramid
R=(1-t).*r;
S=(1-t).*s;
T=t;
C=(1-t).^2.*Cr.*Cs.*Ct;


%% demo
% % Get nodes on a unit hexahedron
% n=5;
% [r, Cr]=GaussLobattoR(n, -1, 1);
% [s, Cs]=GaussLobattoR(n+1, -1, 1);
% [t, Ct]=GaussLobattoR(n+2, 0, 1);
% [R, S, T, C]=meshgridpyramc({r, s, t}, {Cr, Cs, Ct});
% 
% % Symbbolic expression
% syms r s t
% Eq=1+r^6+s^3+t^2;
% fh=@(r, s, t) 1+ r.^6+s.^3+t.^2;
% 
% iEq=int(Eq, r, t-1, 1-t);
% iEq=int(iEq, s, t-1, 1-t);
% iEq=double(int(iEq, t, 0, 1));
% 
% % Numerical intergation
% Z=fh(R, S, T);
% iZ=sum(C(:).*Z(:));
% 
% % Draw the domain
% hexasurfg(R, S, T);
end







