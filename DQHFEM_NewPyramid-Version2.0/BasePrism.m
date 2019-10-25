function [G,Gr,Gs,Gt]=BasePrism(r, s, t, number)

%% 求取三棱柱基的函数值和三个偏导数值
% number为结构体：储存三棱柱单元边、面、体上的阶次
% 可缺省调用[G]=BasePrism(r,s,t,number )
%
NE=number.edge;FA=number.face;n=number.H;
ne=[NE(1:3)+2;NE(7:9)+2;NE(4:6)];
nft=[FA{1},FA{5}];
nfq=[FA{2};FA{3};FA{4}];
% Numbers of nodes on faces
nb1=[ne(1,1)-2, ne(2,1)-2, ne(3,1), ne(3,2)];
nb2=[ne(1,2)-2, ne(2,2)-2, ne(3,3), ne(3,2)];
nb3=[ne(1,3)-2, ne(2,3)-2, ne(3,1), ne(3,3)];
nbf(1)=sum(ne(1,:))-3; tnfs(1)=nbf(1)+nft(1)*(nft(1)+1)/2;
nbf(2)=sum(ne(2,:))-3; tnfs(2)=nbf(2)+nft(2)*(nft(2)+1)/2;
nbf(3) = sum(nb1)+4; tnfs(3)= nbf(3)+nfq(1,1)*nfq(1,2);
nbf(4) = sum(nb2)+4; tnfs(4)= nbf(4)+nfq(2,1)*nfq(2,2);
nbf(5) = sum(nb3)+4; tnfs(5)= nbf(5)+nfq(3,1)*nfq(3,2);
tnf=sum(sum(ne))-6+nft(1)*(nft(1)+1)/2+nft(2)*(nft(2)+1)/2+...
    nfq(1,1)*nfq(1,2)+nfq(2,1)*nfq(2,2)+nfq(3,1)*nfq(3,2);
tn=tnf+n(2)*n(1)*(n(1)+1)/2;

% Prepare matrices
G=zeros(length(r(:)), tn);
if nargout>1
    Gx=G; Gy=G; Gz=G;
end

% Local basis on the triangles, quadrileterals and interior
if nargout==1
    G1=TrigBasisC0(r, s, ne(1,:), nft(1));
    G2=TrigBasisC0(r, s, ne(2,:), nft(2));
    E1=TrigBasisEi(r, s, nfq(1,1), [2, 3, 2, 0], [1,1,0]);
    E2=TrigBasisEp(r, s, nfq(2,1), [0, 5, 2, 2], [0,1,1]);
    E3=TrigBasisEi(r, s, nfq(3,1), [2, 3, 0, 2], [1,0,1]);
    Gt=TrigBasisC0(r, s, [2,2,2], n(1));
    Gq=HierBasisJC0(max([ne(3,:), nfq(:,2)', n(2)]), t(:), 0, 1);
elseif nargout>1
    [G1, jac1]=TrigBasisC0(r, s, ne(1,:), nft(1));
    [G2, jac2]=TrigBasisC0(r, s, ne(2,:), nft(2));
    [E1, jact1]=TrigBasisEi(r, s, nfq(1,1), [2, 3, 2, 0], [1,1,0]);
    [E2, jact2]=TrigBasisEp(r, s, nfq(2,1), [0, 5, 2, 2], [0,1,1]);
    [E3, jact3]=TrigBasisEi(r, s, nfq(3,1), [2, 3, 0, 2], [1,0,1]);
    [Gt, jact]=TrigBasisC0(r, s, [2,2,2], n(1));
    [Gq, jacq]=HierBasisJC0(max([ne(3,:), nfq(:,2)', n(2)]), t(:), 0, 1);
end

% Global basis on the 6 cornors
G(:,1:6)=[(1-r(:)-s(:)).*(1-t(:)), r(:).*(1-t(:)), s(:).*(1-t(:)), (1-r(:)-s(:)).*t(:), r(:).*t(:), s(:).*t(:)];
if nargout>1
    Gx(:,1:6)=[-(1-t(:)), 1-t(:), 0*t(:), -t(:), t(:), 0*t(:)];
    Gy(:,1:6)=[-(1-t(:)), 0*t(:), 1-t(:), -t(:), 0*t(:),t(:)];
    Gz(:,1:6)=[-(1-r(:)-s(:)), -r(:), -s(:), (1-r(:)-s(:)), r(:), s(:)];
end

% Global basis on edge 1
p=6;
for i=1:ne(1,1)-2
    G(:,p+i)=G1(:,1+i).*(1-t(:)); 
    if nargout>1
        Gx(:,p+i)=jac1{1}(:,1+i).*(1-t(:)); 
        Gy(:,p+i)=jac1{2}(:,1+i).*(1-t(:)); 
        Gz(:,p+i)=-G1(:,1+i); 
    end
end
p=p+ne(1,1)-2;

% Global basis on edge 2
for i=1:ne(1,2)-2
    G(:,p+i)=G1(:,ne(1,1)+i).*(1-t(:)); 
    if nargout>1
        Gx(:,p+i)=jac1{1}(:,ne(1,1)+i).*(1-t(:)); 
        Gy(:,p+i)=jac1{2}(:,ne(1,1)+i).*(1-t(:)); 
        Gz(:,p+i)=-G1(:,ne(1,1)+i); 
    end
end
p=p+ne(1,2)-2;

% Global basis on edge 3
for i=1:ne(1,3)-2
    G(:,p+i)=G1(:,sum(ne(1,1:2))-1+i).*(1-t(:)); 
    if nargout>1
        Gx(:,p+i)=jac1{1}(:,sum(ne(1,1:2))-1+i).*(1-t(:)); 
        Gy(:,p+i)=jac1{2}(:,sum(ne(1,1:2))-1+i).*(1-t(:)); 
        Gz(:,p+i)=-G1(:,sum(ne(1,1:2))-1+i);
    end
end
p=p+ne(1,3)-2;

% Global basis on edge 4
for i=1:ne(2,1)-2
    G(:,p+i)=G2(:,1+i).*t(:); 
    if nargout>1
        Gx(:,p+i)=jac2{1}(:,1+i).*t(:); 
        Gy(:,p+i)=jac2{2}(:,1+i).*t(:); 
        Gz(:,p+i)=G2(:,1+i); 
    end
end
p=p+ne(2,1)-2;

% Global basis on edge 5
for i=1:ne(2,2)-2
    G(:,p+i)=G2(:,ne(2,1)+i).*t(:); 
    if nargout>1
        Gx(:,p+i)=jac2{1}(:,ne(2,1)+i).*t(:); 
        Gy(:,p+i)=jac2{2}(:,ne(2,1)+i).*t(:); 
        Gz(:,p+i)=G2(:,ne(2,1)+i); 
    end
end
p=p+ne(2,2)-2;

% Global basis on edge 6
for i=1:ne(2,3)-2
    G(:,p+i)=G2(:,sum(ne(2,1:2))-1+i).*t(:); 
    if nargout>1
        Gx(:,p+i)=jac2{1}(:,sum(ne(2,1:2))-1+i).*t(:); 
        Gy(:,p+i)=jac2{2}(:,sum(ne(2,1:2))-1+i).*t(:); 
        Gz(:,p+i)=G2(:,sum(ne(2,1:2))-1+i); 
    end
end
p=p+ne(2,3)-2;

% Global basis on edge 7
for i=1:ne(3,1)
    G(:,p+i)=Gq(:,2+i).*(1-r(:)-s(:)); 
    if nargout>1
        Gx(:,p+i)=-Gq(:,2+i); 
        Gy(:,p+i)=-Gq(:,2+i); 
        Gz(:,p+i)=jacq(:,2+i).*(1-r(:)-s(:)); 
    end
end
p=p+ne(3,1);

% Global basis on edge 8
for i=1:ne(3,2)
    G(:,p+i)=Gq(:,2+i).*r(:); 
    if nargout>1
        Gx(:,p+i)=Gq(:,2+i); 
        Gy(:,p+i)=0*Gq(:,2+i); 
        Gz(:,p+i)=jacq(:,2+i).*r(:); 
    end
end
p=p+ne(3,2);

% Global basis on edge 9
for i=1:ne(3,3)
    G(:,p+i)=Gq(:,2+i).*s(:); 
    if nargout>1
        Gx(:,p+i)=0*Gq(:,2+i); 
        Gy(:,p+i)=Gq(:,2+i); 
        Gz(:,p+i)=jacq(:,2+i).*s(:); 
    end
end
p=p+ne(3,3);

% Global basis on face 1
for i=1:tnfs(1)-nbf(1)
    G(:,p+i)=G1(:,nbf(1)+i).*(1-t(:)); 
    if nargout>1
        Gx(:,p+i)=jac1{1}(:,nbf(1)+i).*(1-t(:)); 
        Gy(:,p+i)=jac1{2}(:,nbf(1)+i).*(1-t(:)); 
        Gz(:,p+i)=-G1(:,nbf(1)+i); 
    end
end
p=p+tnfs(1)-nbf(1);

% Global basis on face 2
for i=1:tnfs(2)-nbf(2)
    G(:,p+i)=G2(:,nbf(2)+i).*t(:); 
    if nargout>1
        Gx(:,p+i)=jac2{1}(:,nbf(2)+i).*t(:); 
        Gy(:,p+i)=jac2{2}(:,nbf(2)+i).*t(:); 
        Gz(:,p+i)=G2(:,nbf(2)+i); 
    end
end
p=p+tnfs(2)-nbf(2);

% Global basis on face 3
for j=1:nfq(1,2) 
    for i=1:nfq(1,1) 
        G(:,p+(j-1)*nfq(1,1)+i)=E1(:,i).*Gq(:,j+2); 
        if nargout>1
            Gx(:,p+(j-1)*nfq(1,1)+i)=jact1{1}(:,i).*Gq(:,j+2); 
            Gy(:,p+(j-1)*nfq(1,1)+i)=jact1{2}(:,i).*Gq(:,j+2); 
            Gz(:,p+(j-1)*nfq(1,1)+i)=E1(:,i).*jacq(:,j+2); 
        end
    end
end
p=p+tnfs(3)-nbf(3);

% Global basis on face 4
for j=1:nfq(2,2) 
    for i=1:nfq(2,1) 
        G(:,p+(j-1)*nfq(2,1)+i)=E2(:,i).*Gq(:,j+2); 
        if nargout>1
            Gx(:,p+(j-1)*nfq(2,1)+i)=jact2{1}(:,i).*Gq(:,j+2); 
            Gy(:,p+(j-1)*nfq(2,1)+i)=jact2{2}(:,i).*Gq(:,j+2); 
            Gz(:,p+(j-1)*nfq(2,1)+i)=E2(:,i).*jacq(:,j+2); 
        end
    end
end
p=p+tnfs(4)-nbf(4);

% Global basis on face 5
for j=1:nfq(3,2) 
    for i=1:nfq(3,1) 
        G(:,p+(j-1)*nfq(3,1)+i)=E3(:,i).*Gq(:,j+2); 
        if nargout>1
            Gx(:,p+(j-1)*nfq(3,1)+i)=jact3{1}(:,i).*Gq(:,j+2); 
            Gy(:,p+(j-1)*nfq(3,1)+i)=jact3{2}(:,i).*Gq(:,j+2); 
            Gz(:,p+(j-1)*nfq(3,1)+i)=E3(:,i).*jacq(:,j+2); 
        end
    end
end
p=p+tnfs(5)-nbf(5);

% Global basis inside the element
nt=n(1)*(n(1)+1)/2;
for j=1:n(2)
    for i=1:nt
        p=p+1;
        G(:,p)=Gt(:,3+i).*Gq(:,j+2); 
        if nargout>1
            Gx(:,p)=jact{1}(:,3+i).*Gq(:,j+2); 
            Gy(:,p)=jact{2}(:,3+i).*Gq(:,j+2); 
            Gz(:,p)=Gt(:,3+i).*jacq(:,j+2);      
        end
    end
end

if nargout>1
    jac={Gx, Gy, Gz};
    Gr=jac{1};Gs=jac{2};Gt=jac{3};
end

%% demo - draw basis
% % The number of nodes on triangular prism edges (ne), 
% % faces (nft, nfq), and the number inside hierarchical nodes (n)
% % The number of nodes for plotting (M, N)
% ne=[3,3,3; 3,3,3; 2,2,2]; nft=[2,2]; nfq=[2,2; 2,2; 2,2]; n=[2, 2];
% M=12; N=M;
% 
% % Get the interpolation nodes on a unit triangular prism
% [x, y, z, tri]=prismnodes(ne, nft, nfq, n);
% 
% % Get basis functions on a unit triangular prism and their derivatives
% [r, s, t]=meshgridprismc(M, N);
% [G, jac]=PrismBasis(r, s, t, ne, nft, nfq, n);
% 
% % Plot the results
% j=24; 
% Gj=reshape(G(:,j), M, M, N);
% figure; hold on;
% prismsurfg(r, s, t, Gj);
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% shading interp;
% colormap jet;
% colorbar; alpha(0.2);
% axis equal;
% view(3);


%% demo - derivatives of basis
% % The number of nodes on triangular prism edges (ne), 
% % faces (nft, nfq), and the number inside hierarchical nodes (n)
% % The number of nodes for plotting (M, N)
% ne=[3,3,3; 3,3,3; 2,2,2]; nft=[2,2]; nfq=[2,2; 2,2; 2,2]; n=[2, 2];
% M=12; N=M;
% 
% % The number of points covered by RBF
% nup=4*N; % this parameter can be optimized ...
% 
% % Get the interpolation nodes on a unit triangular prism
% [x, y, z]=prismnodes(ne, nft, nfq, n);
% 
% % Get basis functions on a unit triangular prism and their derivatives
% [r, s, t, tri]=prismnodes(M*ones(3), M*ones(1,2), M*ones(3,2), [N,N]);
% [G, jac]=PrismBasis(r, s, t, ne, nft, nfq, n);
% Gx=jac{1}; Gy=jac{2}; Gz=jac{3}; 
% 
% % Radial basis 3D
% pnts=[r(:)'; s(:)'; t(:)']; 
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% Gxi=Dx*G;
% Gyi=Dy*G;
% Gzi=Dz*G;
% 
% % Plot the results
% j=42;
% figure; hold on;
% for k=1:5+N
%     trisurf(tri{k}, r, s, t, G(:,j));
% end
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% colormap jet;
% shading interp;
% colorbar; alpha(0.3);
% axis equal;
% view(3);
% 
% figure; hold on;
% for k=1:5+N
%     trisurf(tri{k}, r, s, t, Gx(:,j));
% end
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% colormap jet;
% shading interp;
% colorbar; alpha(0.3);
% axis equal;
% view(3);
% title('Exact');
% 
% figure; hold on;
% for k=1:5+N
%     trisurf(tri{k}, r, s, t, Gxi(:,j));
% end
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% colormap jet;
% shading interp;
% colorbar; alpha(0.3);
% axis equal;
% view(3);
% title('RBF');


%% demo - interpolation
% % The number of nodes on triangular prism edges (ne), 
% % faces (nft, nfq), and the number inside hierarchical nodes (n)
% % The number of nodes for plotting (M, N)
% K=3; ne=[3,3,3; 3,3,3; 1,1,1]+K; nft=[2,2]+K; 
% nfq=[2,2; 2,2; 2,2]+K; n=[2, 2]+K;
% M=15; N=M;
% 
% % Get the interpolation nodes on a unit triangular prism
% [x, y, z]=prismnodes(ne, nft, nfq, n);
% 
% % Get basis functions on a unit triangular prism and their derivatives
% Q=PrismBasis(x, y, z, ne, nft, nfq, n); iQ=inv(Q);
% [r, s, t, tri]=prismnodes(M*ones(3), M*ones(1,2), M*ones(3,2), [N,N]);
% [G, jac]=PrismBasis(r, s, t, ne, nft, nfq, n);
% G=G*iQ; Gx=jac{1}*iQ; Gy=jac{2}*iQ; Gz=jac{3}*iQ; 
% 
% % Function values and derivatives
% fun=@(x, y, z) cos(x).*cos(y).*cos(z);
% funx=@(x, y, z) -sin(x).*cos(y).*cos(z);
% funy=@(x, y, z) -cos(x).*sin(y).*cos(z);
% funz=@(x, y, z) -cos(x).*cos(y).*sin(z);
% F0=fun(x, y, z); F=fun(r, s, t); Fx=funx(r, s, t); 
% Fy=funy(r, s, t); Fz=funz(r, s, t); 
% Fi=G*F0; Fxi=Gx*F0; Fyi=Gy*F0; Fzi=Gz*F0; 
% 
% % Plot the results
% figure; hold on;
% for k=1:5+N
%     trisurf(tri{k}, r, s, t, Fz);
% end
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% colormap jet;
% shading interp;
% colorbar; alpha(0.3);
% axis equal;
% view(3);
% 
% figure; hold on;
% for k=1:5+N
%     trisurf(tri{k}, r, s, t, Fz-Fzi);
% end
% plot3(x, y, z, 'r.', 'MarkerSize', 20);
% colormap jet;
% shading interp;
% colorbar; alpha(0.3);
% axis equal;
% view(3);
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G, jac, hess]=TrigBasisC0(x, y, mpn, K)

% TrigBasisF : Get C0 basis functions and their derivatives on a unit triangle
% 
% Calling Sequences:
%
%     G=TrigBasisC0(x, y, mpn, K)
%
%     [G, jac]=TrigBasisC0(x, y, mpn, K)
%
%     [G, jac, hess]=TrigBasisC0(x, y, mpn, K)
% 
% INPUTS:
% 
%     x, y   -    Vectors of parametric coordinates of a unit triangle.
%
%     mpn = [M, P, N] - the number of nodes on triangle edges
%
%     K - The number of basis on the unit triangle.
%
% OUTPUT:
% 
%      G   -  Basis functions
%
%      jac - First derivatives of the Face functions
%
%      hess - Second derivatives of the Face functions
% 

% Basis on the boundary of a unit triangle
M=mpn(1); P=mpn(2); N=mpn(3); 
x=x(:); y=y(:);
G=zeros(length(x), M+N+P-3); 

if nargout==1
    % Conor basis
    G(:, [1,M,M+P-1])=[1-x-y, x, y];
    
    % Basis on edge 1
    pp=2:M-1;
    G(:, pp)=TrigBasisEi(x, y, M-2, [2, 3, 2, 0], [1,1,0]);
    
    % Basis on edge 2
    pp=M+1:M+P-2;
    G(:, pp)=TrigBasisEp(x, y, P-2, [0, 5, 2, 2], [0,1,1]);
    
    % Basis on edge 3
    pp=M+P:M+P+N-3;
    G(:, pp)=TrigBasisEi(x, y, N-2, [2, 3, 0, 2], [1,0,1]);
    
    % Basis on the face of the element
    F=TrigBasisF(x, y, K, [2, 5, 2, 2], [1,1,1]);
    
    % All basis (edges and face basis)
    G=[G, F]; 
elseif nargout==2
    Gx=G; Gy=G;
    
    % Conor basis
    wn=ones(size(x));
    G(:, [1,M,M+P-1])=[1-x-y, x, y];
    Gx(:, [1,M,M+P-1])=[-wn, wn, 0*wn];
    Gy(:, [1,M,M+P-1])=[-wn, 0*wn, wn];
    
    % Basis on edge 1
    pp=2:M-1;
    [E1, jac1]=TrigBasisEi(x, y, M-2, [2, 3, 2, 0], [1,1,0]);
    G(:, pp)=E1; Gx(:, pp)=jac1{1}; Gy(:, pp)=jac1{2}; 
    
    % Basis on edge 2
    pp=M+1:M+P-2;
    [E2, jac2]=TrigBasisEp(x, y, P-2, [0, 5, 2, 2], [0,1,1]);
    G(:, pp)=E2; Gx(:, pp)=jac2{1}; Gy(:, pp)=jac2{2}; 
    
    % Basis on edge 3
    pp=M+P:M+P+N-3;
    [E3, jac3]=TrigBasisEi(x, y, N-2, [2, 3, 0, 2], [1,0,1]);
    G(:, pp)=E3; Gx(:, pp)=jac3{1}; Gy(:, pp)=jac3{2}; 
    
    % Basis on the face of the element
    [F, jac]=TrigBasisF(x, y, K, [2, 5, 2, 2], [1,1,1]);
    
    % All basis (edges and face basis)
    G=[G, F]; jac{1}=[Gx, jac{1}]; jac{2}=[Gy, jac{2}]; 
elseif nargout==3
    Gx=G; Gy=G; Gxx=G; Gxy=G; Gyy=G;
    
    % Conor basis
    wn=ones(size(x));
    G(:, [1,M,M+P-1])=[1-x-y, x, y];
    Gx(:, [1,M,M+P-1])=[-wn, wn, 0*wn];
    Gy(:, [1,M,M+P-1])=[-wn, 0*wn, wn];
    
    % Basis on edge 1
    pp=2:M-1;
    [E1, jac1, hess1]=TrigBasisEi(x, y, M-2, [2, 3, 2, 0], [1,1,0]);
    G(:, pp)=E1; Gx(:, pp)=jac1{1}; Gy(:, pp)=jac1{2}; 
    Gxx(:, pp)=hess1{1,1}; Gxy(:, pp)=hess1{1,2}; 
    Gyy(:, pp)=hess1{2,2}; 
    
    % Basis on edge 2
    pp=M+1:M+P-2;
    [E2, jac2, hess2]=TrigBasisEp(x, y, P-2, [0, 5, 2, 2], [0,1,1]);
    G(:, pp)=E2; Gx(:, pp)=jac2{1}; Gy(:, pp)=jac2{2}; 
    Gxx(:, pp)=hess2{1,1}; Gxy(:, pp)=hess2{1,2}; 
    Gyy(:, pp)=hess2{2,2}; 
    
    % Basis on edge 3
    pp=M+P:M+P+N-3;
    [E3, jac3, hess3]=TrigBasisEi(x, y, N-2, [2, 3, 0, 2], [1,0,1]);
    G(:, pp)=E3; Gx(:, pp)=jac3{1}; Gy(:, pp)=jac3{2}; 
    Gxx(:, pp)=hess3{1,1}; Gxy(:, pp)=hess3{1,2}; 
    Gyy(:, pp)=hess3{2,2}; 
    
    % Basis on the face of the element
    [F, jac, hess]=TrigBasisF(x, y, K, [2, 5, 2, 2], [1,1,1]);
    
    % All basis (edges and face basis)
    G=[G, F]; jac{1}=[Gx, jac{1}]; jac{2}=[Gy, jac{2}]; 
    hess{1,1}=[Gxx, hess{1,1}];
    hess{1,2}=[Gxy, hess{1,2}];
    hess{2,2}=[Gyy, hess{2,2}];
    hess{2,1}=hess{1,2};
end

%% demo - Test for basis
% % Number of basis on triangle edges (M, P, N) and face (Ns)
% M=6; P=5; N=4; Ns=4; n=20;
% 
% % Numbers of hierarchical basis
% NB=M+N+P-3; TN=NB+Ns*(Ns+1)/2;
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% pnts=[x(:)'; y(:)']; 
% 
% % Get the basis on triangle
% [F, jac, hess]=TrigBasisC0(x, y, [M, P, N], Ns);
% Fx=jac{1}; Fy=jac{2}; Fxx=hess{1,1};
% Fxy=hess{1,2}; Fyy=hess{2,2};
% 
% % Interpolate the first derivetives by RBF
% [G, Dx, Dy]=rbfbasis2d(pnts, pnts, 3*n);
% Fxi=Dx*F;
% Fyi=Dy*F;
% Fxxi=Dx*Fx;
% Fxyi=Dy*Fx;
% Fyyi=Dy*Fy;
% 
% k=18;
% figure; hold on;
% trisurf(tri,x,y,Fxy(:,k)); 
% shading interp;
% view(3);
% title('Direct');
% 
% figure; hold on;
% trisurf(tri,x,y,Fxyi(:,k)); 
% shading interp;
% view(3);
% title('RBF');

%% demo - Test for interpolation in natrual coordinates 1
% % The number of basis on triangle edges (M, P, N) and face (Ns)
% %  and the number of nodes to plot the basis (n)
% Ns=11; M=Ns+2; P=M+1; N=M+2; 
% 
% % Numbers of hierarchical basis
% NB=M+N+P-3; tn=NB+Ns*(Ns+1)/2;
% 
% % Get the integration nodes and weights on a unit triangle
% [x, y, c]=GaussTrigH([M, P, N], Ns);
% tri=tridelaunayh(x, y, [M, P, N], Ns);
% 
% % Interpolation
% [F, Fx, Fy, Fxx, Fxy, Fyy]=funz(x, y);
% [G, jac, hess]=TrigBasisC0(x, y, [M, P, N], Ns); 
% Q=inv(G); 
% Gx=jac{1}*Q; 
% Gy=jac{2}*Q;
% Gxx=hess{1,1}*Q; 
% Gxy=hess{1,2}*Q; 
% Gyy=hess{2,2}*Q; 
% Fxi=Gx*F;
% Fyi=Gy*F;
% Fxxi=Gxx*F;
% Fxyi=Gxy*F;
% Fyyi=Gyy*F;
% 
% % Plot the results
% figure; hold on;
% trisurf(tri, x, y, Fyy); 
% view(3); title('Fyy');
% 
% figure; hold on;
% trisurf(tri, x, y, Fyy-Fyyi); 
% view(3); title('Fyyi');


%% demo - Test for interpolation in natrual coordinates 2
% % The number of basis on triangle edges (M, P, N) and face (Ns)
% %  and the number of nodes to plot the basis (n)
% % Added new node numbers (k)
% Ns=11; M=Ns+2; P=M+1; N=M+2; k=2;
% 
% % Numbers of hierarchical basis
% NB=M+N+P-3; tn=NB+Ns*(Ns+1)/2;
% 
% % Get the integration nodes and weights on a unit triangle
% [s, t]=GaussTrigH([M, P, N], Ns);
% [x, y]=GaussTrigH([M, P, N]+k, Ns+k);
% tri=tridelaunayh(x, y, [M, P, N]+k, Ns+k);
% 
% % Interpolation
% F0=funf(s, t);
% [F, Fx, Fy, Fxx, Fxy, Fyy]=funf(x, y);
% Q=TrigBasisC0(s, t, [M, P, N], Ns); 
% [G, jac, hess]=TrigBasisC0(x, y, [M, P, N], Ns); 
% Q=inv(Q); 
% G=G*Q; 
% Gx=jac{1}*Q; 
% Gy=jac{2}*Q;
% Gxx=hess{1,1}*Q; 
% Gxy=hess{1,2}*Q; 
% Gyy=hess{2,2}*Q; 
% Fi=G*F0;
% Fxi=Gx*F0;
% Fyi=Gy*F0;
% Fxxi=Gxx*F0;
% Fxyi=Gxy*F0;
% Fyyi=Gyy*F0;
% 
% % Plot the results
% figure; hold on;
% trisurf(tri, x, y, F); 
% view(3); title('F');
% 
% figure; hold on;
% trisurf(tri, x, y, F-Fi); 
% view(3); title('Fi');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, jac, hess, num]=TrigBasisEi(x, y, N, jc, pe)

% TrigBasisEp : Get edge functions i on a unit triangle and their derivatives
% 
% Calling Sequences:
%
%     F=TrigBasisEi(x, y, N, [a,b,c,d], pe)
%
%     [F, jac]=TrigBasisEi(x, y, N, [a,b,c,d], pe)
%
%     [F, jac, hess, num]=TrigBasisEi(x, y, N, [a,b,c,d], pe)
% 
% INPUTS:
% 
%     x, y   -    Vectors of parametric coordinates of a unit triangle
%
%     N - The number of basis on the unit triangle
%
%     jc=[a,b,c,d]  -  Parameters of Jacobi polynomials on the unit triangle.
%                     The parameters c and d are not needed in this case.
%
%     pe  -  Coefficients of the power of the polynomial of the unit triangle
%
% OUTPUT:
% 
%      F   -  Face functions
%
%      jac - First derivatives of the Face functions
%
%      hess - Second derivatives of the Face functions
%
%      num - A matrix of indexes of the face functions
% 

% Prepare variables and matrices
x=x(:); y=y(:);
a=jc(1); b=jc(2); 
F=zeros(length(x), N); 
num=zeros(N, 2);
if nargout>=2
    Fx=F; Fy=F; 
end
if nargout>=3
    Fxx=F; Fxy=F; Fyy=F;
end

% Get the basis on triangle
if nargout==1
    Pt=polytrider(x, y, pe(1), pe(2), pe(3));
    dJ=JacobiDerCellA(N-1, a, b, 2*(x+y) - 1, 0);
    Pi=cell2mat(dJ(:,1)'); 
elseif nargout==2
    [Pt, jact]=polytrider(x, y, pe(1), pe(2), pe(3));   
    dJ=JacobiDerCellA(N-1, a, b, 2*(x+y) - 1, 1);
    Pi=cell2mat(dJ(:,1)'); 
    dPi=cell2mat(dJ(:,2)'); 
elseif nargout>=3
    [Pt, jact, hesst]=polytrider(x, y, pe(1), pe(2), pe(3));
    dJ=JacobiDerCellA(N-1, a, b, 2*(x+y) - 1, 2);
    Pi=cell2mat(dJ(:,1)'); 
    dPi=cell2mat(dJ(:,2)'); 
    ddPi=cell2mat(dJ(:,3)'); 
end
t=1; m=1;
for n=1:N-m+1
    if nargout==2
        jach={2*dPi(:,n), 2*dPi(:,n)};
    elseif nargout>=3
        jach={2*dPi(:,n), 2*dPi(:,n)};
        hessh={4*ddPi(:,n), 4*ddPi(:,n); 4*ddPi(:,n), 4*ddPi(:,n)};
    end

    % The value of the basis
    Ph=Pi(:,n);
    F(:,t)=Pt(:).*Ph;

    % First derivative of the basis
    if nargout>=2
        Fx(:,t)=Pt(:).*jach{1}+Ph.*jact{1}(:);
        Fy(:,t)=Pt(:).*jach{2}+Ph.*jact{2}(:);
    end

    % Second derivative of the basis
    if nargout>=3
        Fxx(:,t)=Pt(:).*hessh{1,1}+2*jact{1}(:).*jach{1}+Ph.*hesst{1,1}(:);
        Fyy(:,t)=Pt(:).*hessh{2,2}+2*jact{2}(:).*jach{2}+Ph.*hesst{2,2}(:);
        Fxy(:,t)=jact{1}(:).*jach{2}+Pt(:).*hessh{1,2}+Ph.*hesst{1,2}(:)+jach{1}.*jact{2}(:);
    end
    num(t,1)=m; num(t,2)=n; 
    t=t+1;
end
if nargout>=2
    jac={Fx, Fy};
end
if nargout>=3
    hess={Fxx, Fxy; Fxy, Fyy};
end
%% demo
% % Node numbers (n) and the number of basis (N>=0)
% % Parameters of Jacobi polynomials (a, b, c, d)
% n=17; N=5; a=2; b=3; c=2; d=0; pe=[1,0,1];
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% pnts=[x(:)'; y(:)']; 
% 
% % Get the basis on triangle
% [F, jac, hess, num]=TrigBasisEi(x, y, N, [a, b, c, d], pe);
% Fx=jac{1}; Fy=jac{2}; Fxx=hess{1,1};
% Fxy=hess{1,2}; Fyy=hess{2,2};
% 
% % Interpolate the first derivetives by RBF
% [G, Dx, Dy]=rbfbasis2d(pnts, pnts, 3*n);
% Fxi=Dx*F;
% Fyi=Dy*F;
% Fxxi=Dx*Fx;
% Fxyi=Dy*Fx;
% Fyyi=Dy*Fy;
% 
% k=3;
% figure; hold on;
% trisurf(tri,x,y,Fy(:,k)); 
% shading interp;
% view(3);
% title('Direct');
% 
% figure; hold on;
% trisurf(tri,x,y,Fyi(:,k)); 
% shading interp;
% view(3);
% title('RBF');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, jac, hess, num]=TrigBasisEp(x, y, N, jc, pe)

% TrigBasisF : Get edge functions p on a unit triangle and their derivatives
% 
% Calling Sequences:
%
%     F=TrigBasisEp(x, y, N, [a,b,c,d], pe)
%
%     [F, jac]=TrigBasisEp(x, y, N, [a,b,c,d], pe)
%
%     [F, jac, hess, num]=TrigBasisEp(x, y, N, [a,b,c,d], pe)
% 
% INPUTS:
% 
%     x, y   -    Vectors of parametric coordinates of a unit triangle
%
%     N - The number of basis on the unit triangle
%
%     jc=[a,b,c,d]  -  Parameters of Jacobi polynomials on the unit triangle.
%                     The parameters a and b are not needed in this case.
%
%     pe  -  Coefficients of the power of the polynomial of the unit triangle
%
% OUTPUT:
% 
%      F   -  Face functions
%
%      jac - First derivatives of the Face functions
%
%      hess - Second derivatives of the Face functions
%
%      num - A matrix of indexes of the face functions
% 

% Prepare variables and matrices
x=x(:); y=y(:);
c=jc(3); d=jc(4); 
F=zeros(length(x), N); 
num=zeros(N, 2);
if nargout>=2
    Fx=F; Fy=F; 
end
if nargout>=3
    Fxx=F; Fxy=F; Fyy=F;
end

% Get the basis on triangle
if nargout==1
    Ph=JacobiTrig(N-1, c, d, x, y);
    Pt=polytrider(x, y, pe(1), pe(2), pe(3));    
elseif nargout==2
    [Ph, jach]=JacobiTrig(N-1, c, d, x, y);
    [Pt, jact]=polytrider(x, y, pe(1), pe(2), pe(3));   
elseif nargout>=3
    [Ph, jach, hessh]=JacobiTrig(N-1, c, d, x, y);
    [Pt, jact, hesst]=polytrider(x, y, pe(1), pe(2), pe(3));
end
t=1; 
for m=1:N    
    % The value of the basis
    F(:,t)=Pt(:).*Ph(:,m);

    % First derivative of the basis
    if nargout>=2
        Fx(:,t)=Pt(:).*jach{1}(:,m)+Ph(:,m).*jact{1}(:);
        Fy(:,t)=Pt(:).*jach{2}(:,m)+Ph(:,m).*jact{2}(:);
    end

    % Second derivative of the basis
    if nargout>=3
        Fxx(:,t)=Pt(:).*hessh{1,1}(:,m)+2*jact{1}(:).*jach{1}(:,m)+Ph(:,m).*hesst{1,1}(:);
        Fyy(:,t)=Pt(:).*hessh{2,2}(:,m)+2*jact{2}(:).*jach{2}(:,m)+Ph(:,m).*hesst{2,2}(:);
        Fxy(:,t)=jact{1}(:).*jach{2}(:,m)+Pt(:).*hessh{1,2}(:,m)+Ph(:,m).*hesst{1,2}(:)+jach{1}(:,m).*jact{2}(:);
    end
    num(t,1)=m; num(t,2)=1; 
    t=t+1;
end
if nargout>=2
    jac={Fx, Fy};
end
if nargout>=3
    hess={Fxx, Fxy; Fxy, Fyy};
end

%% demo
% % Node numbers (n) and the number of basis (N>=0)
% % Parameters of Jacobi polynomials (a, b, c, d)
% n=17; N=5; a=0; b=5; c=2; d=2; pe=[0,1,1];
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% pnts=[x(:)'; y(:)']; 
% 
% % Get the basis on triangle
% [F, jac, hess, num]=TrigBasisEp(x, y, N, [a, b, c, d], pe);
% Fx=jac{1}; Fy=jac{2}; Fxx=hess{1,1};
% Fxy=hess{1,2}; Fyy=hess{2,2};
% 
% % Interpolate the first derivetives by RBF
% [G, Dx, Dy]=rbfbasis2d(pnts, pnts, 3*n);
% Fxi=Dx*F;
% Fyi=Dy*F;
% Fxxi=Dx*Fx;
% Fxyi=Dy*Fx;
% Fyyi=Dy*Fy;
% 
% k=3;
% figure; hold on;
% trisurf(tri,x,y,Fx(:,k)); 
% shading interp;
% view(3);
% title('Direct');
% 
% figure; hold on;
% trisurf(tri,x,y,Fxi(:,k)); 
% shading interp;
% view(3);
% title('RBF');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [G, Gx, Gxx, Gxxx]=HierBasisJC0(M, x, a ,b)

% Get 1D hierarchical Jacobi basis for C0 problem
% 
%  Calling sequence:  
% 
%    G=HierBasisJC0(M, x)
%    G=HierBasisJC0(M, x, a ,b) 
%
%    [G, Gx]=HierBasisJC0(M, x)
%    [G, Gx]=HierBasisJC0(M, x, a ,b)
%
%    [G, Gx, Gxx]=HierBasisJC0(M, x)
%    [G, Gx, Gxx]=HierBasisJC0(M, x, a ,b)
%
%    [G, Gx, Gxx, Gxxx]=HierBasisJC0(M, x)
%    [G, Gx, Gxx, Gxxx]=HierBasisJC0(M, x, a ,b)
%   
%  Input :
%
%    M = the number of hierarchical basis
% 
%    x - the integration or interpolation nodes
% 
%    [a, b] - the interval of the domain. Default values are [-1, 1]
% 
%  Output :
% 
%    G, Gx, Gxx, Gxxx - the interpolation matrix and its derivatives
% 

if nargin~=2 && nargin~=4
    error('Error number of inputs');
end
if nargin==4
    x=(x-(b+a)/2)/((b-a)/2);
end

% The 1D basis for C0 problem
Nx=length(x);
G=zeros(Nx, M+2); 
G(:,1)=0.5*(1-x); G(:,2)=0.5*(1+x);
if nargout>1
    Gx=zeros(Nx, M+2);
    Gx(:,1)=-0.5; Gx(:,2)=0.5; 
end
if nargout>2
    Gxx=zeros(Nx, M+2);
end
if nargout>3
    Gxxx=zeros(Nx, M+2);
end
n=nargout-1;
dHJ=HierarJC0(M, x, n);
if M>0
    G(:,3:end)=cell2mat(dHJ(:,1)');
    if nargout>1
        Gx(:,3:end)=cell2mat(dHJ(:,2)');
    end
    if nargout>2
        Gxx(:,3:end)=cell2mat(dHJ(:,3)');
    end
    if nargout>3
        Gxxx(:,3:end)=cell2mat(dHJ(:,4)');
    end
end
if nargin==4
    if nargout>1
        Gx=Gx/((b-a)/2); 
    end
    if nargout>2
        Gxx=Gxx/((b-a)/2)^2; 
    end
    if nargout>3
        Gxxx=Gxxx/((b-a)/2)^3; 
    end
end


%% demo - test of basis
% % Length of the bar (L) and the number of hierarchical basis (M)
% M=20; L=10; 
% 
% % The 1D basis for C0 problem
% a=0; b=L; 
% [x, Cx]=GaussLobattoR(M+4, a, b); 
% [G, Gx, Gxx]=HierBasisJC0(M, x, a ,b);
% 
% % Interpolate the first derivetives by DQM
% Dx=Weighting(x);
% Fxi=Dx*G;
% Fxxi=Dx*Gx;
% 
% % Plot the basis
% k=4;
% figure; plot(x, G(:,k)); title('Basis');
% figure; plot(x, Gx(:,k)); title('First derivative');
% figure; plot(x, Gxx(:,k)); title('Second derivative');


%% demo - vibration of bar
% % Length of the bar (L), Young's modulus (E) and density (rho)
% % The number of hierarchical basis (M)
% M=100; L=10; dM=10; E=71e9; rho=2700;
% 
% % Boundary conditions  (BC): 1 - Clamped, 2 - Free
% BC=[1, 1];
% 
% % The 1D basis for C0 problem
% a=0; b=L; 
% [x, Cx]=GaussLobattoR(M+4, a, b); 
% [G, Gx]=HierBasisJC0(M, x, a ,b);
% 
% % Stiffness and mass matrices
% C=diag(Cx);
% Ke=E*Gx'*C*Gx; Me=rho*G'*C*G;
% 
% % Apply boundary conditions
% bc=true(M+2, 1);
% if BC(1)==1
%     bc(1)=false;
% end
% if BC(2)==1
%     bc(2)=false;
% end
% Ke=Ke(bc, bc); 
% Me=Me(bc, bc);
% 
% % Solve Eigenvalues (SE): 1 - all, 2,3- first a few
% Ns=50; SE=1; 
% [d, V, TN]=SolveEig(Ke, Me, Ns, SE); 
% [sd, I]=sort(d); 
% omg=real(sqrt(sd));
% sdd=omg*L/(pi*sqrt(E/rho));
% 
% % Plot mode
% k=1;
% Vt=zeros(M+2, 1);
% Vt(bc)=V(:, I(k));
% Vt=G*Vt;
% figure; plot(x, Vt);
% 
% % Plot error
% esdd=(1:TN)';
% figure; plot((1:TN)/TN,(sdd(:)-esdd)./esdd)
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F, jac, hess, num]=TrigBasisF(x, y, N, jc, pe)

% TrigBasisF : Get face functions on a unit triangle and their derivatives
% 
% Calling Sequences:
%
%     F=TrigBasisF(x, y, N, [a,b,c,d], pe)
%
%     [F, jac]=TrigBasisF(x, y, N, [a,b,c,d], pe)
%
%     [F, jac, hess, num]=TrigBasisF(x, y, N, [a,b,c,d], pe)
% 
% INPUTS:
% 
%     x, y   -    Vectors of parametric coordinates of a unit triangle
%
%     N - The number of nodes on the unit triangle
%
%     jc=[a,b,c,d]  -  Parameters of Jacobi polynomials on the unit triangle
%
%     pe  -  Coefficients of the power of the polynomial of the unit triangle
%
% OUTPUT:
% 
%      F   -  Face functions
%
%      jac - First derivatives of the Face functions
%
%      hess - Second derivatives of the Face functions
%
%      num - A matrix of indexes of the face functions
% 

% Prepare variables and matrices
a=jc(1); b=jc(2); c=jc(3); d=jc(4); 
F=zeros(length(x), N*(N+1)/2); 
num=zeros(N*(N+1)/2, 2);
if nargout>=2
    Fx=F; Fy=F; 
end
if nargout>=3
    Fxx=F; Fxy=F; Fyy=F;
end

% Get the basis on triangle
if nargout==1
    Pu=JacobiTrig(N-1, c, d, x, y);
    Pt=polytrider(x, y, pe(1), pe(2), pe(3));    
elseif nargout==2
    [Pu, jacu]=JacobiTrig(N-1, c, d, x, y);
    [Pt, jact]=polytrider(x, y, pe(1), pe(2), pe(3));  
elseif nargout>=3
    [Pu, jacu, hessu]=JacobiTrig(N-1, c, d, x, y);
    [Pt, jact, hesst]=polytrider(x, y, pe(1), pe(2), pe(3));
end
t=1;
for m=1:N
    for n=1:N-m+1
        p=m-1; i=n-1; 
        if nargout==1
            Pi=JacobiRecDer(i, a, 2*p+b, 2*(x+y) - 1);
        elseif nargout==2
            [Pi, dPi]=JacobiRecDer(i, a, 2*p+b, 2*(x+y) - 1);
            jaci={2*dPi, 2*dPi};
        elseif nargout>=3
            [Pi, dPi, ddPi]=JacobiRecDer(i, a, 2*p+b, 2*(x+y) - 1);
            jaci={2*dPi, 2*dPi};
            hessi={4*ddPi, 4*ddPi; 4*ddPi, 4*ddPi};
        end
        
        % The value of the basis
        Ph=Pi.*Pu(:,m);
        F(:,t)=Pt(:).*Ph;
        
        % First derivative of the basis
        if nargout>=2
            jach{1}=Pi.*jacu{1}(:,m)+Pu(:,m).*jaci{1};
            jach{2}=Pi.*jacu{2}(:,m)+Pu(:,m).*jaci{2};
            Fx(:,t)=Pt(:).*jach{1}+Ph.*jact{1}(:);
            Fy(:,t)=Pt(:).*jach{2}+Ph.*jact{2}(:);
        end
        
        % Second derivative of the basis
        if nargout>=3
            hessh{1,1}=Pu(:,m).*hessi{1,1}+2*jacu{1}(:,m).*jaci{1}+Pi.*hessu{1,1}(:,m);
            hessh{2,2}=Pu(:,m).*hessi{2,2}+2*jacu{2}(:,m).*jaci{2}+Pi.*hessu{2,2}(:,m);
            hessh{1,2}=jacu{1}(:,m).*jaci{2}+Pu(:,m).*hessi{1,2}+Pi.*hessu{1,2}(:,m)+jacu{2}(:,m).*jaci{1};
            Fxx(:,t)=Pt(:).*hessh{1,1}+2*jact{1}(:).*jach{1}+Ph.*hesst{1,1}(:);
            Fyy(:,t)=Pt(:).*hessh{2,2}+2*jact{2}(:).*jach{2}+Ph.*hesst{2,2}(:);
            Fxy(:,t)=jact{1}(:).*jach{2}+Pt(:).*hessh{1,2}+Ph.*hesst{1,2}(:)+jach{1}.*jact{2}(:);
        end
        num(t,1)=m; num(t,2)=n; 
        t=t+1;
    end
end
if nargout>=2
    jac={Fx, Fy};
end
if nargout>=3
    hess={Fxx, Fxy; Fxy, Fyy};
end


%% demo
% % Node numbers (n) and the number of basis (N>=0)
% % Parameters of Jacobi polynomials (a, b, c, d)
% n=17; N=4; a=2; b=5; c=2; d=2; pe=[1,1,1];
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% pnts=[x(:)'; y(:)']; 
% 
% % Get the basis on triangle
% [F, jac, hess, num]=TrigBasisF(x, y, N, [a, b, c, d], pe);
% Fx=jac{1}; Fy=jac{2}; Fxx=hess{1,1};
% Fxy=hess{1,2}; Fyy=hess{2,2};
% 
% % Interpolate the first derivetives by RBF
% [G, Dx, Dy]=rbfbasis2d(pnts, pnts, 3*n);
% Fxi=Dx*F;
% Fyi=Dy*F;
% Fxxi=Dx*Fx;
% Fxyi=Dy*Fx;
% Fyyi=Dy*Fy;
% 
% k=2;
% figure; hold on;
% trisurf(tri,x,y,Fyy(:,k)); 
% shading interp;
% view(3);
% title('Direct');
% 
% figure; hold on;
% trisurf(tri,x,y,Fyyi(:,k)); 
% shading interp;
% view(3);
% title('RBF');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [p, jac, hess]=polytrider(x, y, a, b, c)

% Get derivatives of polynomials of a unit triangle: (1-x-y)^a*x^b*y^c
% 
% Calling Sequences:
%
%    [p, jac, hess]=polytrider(x, y, a, b, c)
% 
% INPUTS:
% 
%     x, y   -   Coordinates of parametric points.
%
%     a, b, c  - Power of the variables.
%
% OUTPUT:
% 
%     p   -   The value of the expression.
%
%     jac  - Evaluated first derivatives (Jacobian).
% 
%     hess - Evaluated second derivatives (Hessian).
% 

% Get the value of (1-x-y)^a*x^b*y^c
if a<0 || b<0 || c<0
    error('The order of the polynomial cannot be negative!');
end
p=polyeval([1-x-y, x, y], 1, [a,b,c]);

% First derivative
if nargout>=2
    px=polyeval([1-x-y, x, y], -a, [a-1,b,c])+polyeval([1-x-y, x, y], b, [a,b-1,c]);
    py=polyeval([1-x-y, x, y], -a, [a-1,b,c])+polyeval([1-x-y, x, y], c, [a,b,c-1]);
    jac={px, py};
end

% Second derivative
if nargout==3
    pxy=polyeval([1-x-y, x, y], (a-1)*a, [a-2,b,c])+polyeval([1-x-y, x, y], -a*b, [a-1,b-1,c]) ...
            +polyeval([1-x-y, x, y], -a*c, [a-1,b,c-1])+polyeval([1-x-y, x, y], b*c, [a,b-1,c-1]);        
    pxx=polyeval([1-x-y, x, y], (a-1)*a, [a-2,b,c])+2*polyeval([1-x-y, x, y], -a*b, [a-1,b-1,c]) ...
            +polyeval([1-x-y, x, y], (b-1)*b, [a,b-2,c]);
    pyy=polyeval([1-x-y, x, y], (a-1)*a, [a-2,b,c])+2*polyeval([1-x-y, x, y], -a*c, [a-1,b,c-1]) ...
            +polyeval([1-x-y, x, y], (c-1)*c, [a,b,c-2]);
    hess={pxx, pxy; pxy, pyy};
end


%% demo
% % Node numbers (n) and the order of the polynomials (a, b, c)
% n=17; N=6; a=2; b=2; c=2; 
% 
% % The number of points covered by RBF
% nup=3*n; % this parameter can be optimized ...
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% 
% % Get value and derivatives of the polynomial 
% [p, jac, hess]=polytrider(x, y, a, b, c);
% px=jac{1}; py=jac{2};
% pxx=hess{1,1}; pxy=hess{1,2}; pyy=hess{2,2};
% 
% % Plot results
% figure; hold on;
% trisurf(tri,x,y,pyy); 
% view(3);
% title('Results');
% 
% syms s t
% fh=polydeval((1-s-t)^a*s^b*t^c, [s, t], [0,2]);
% figure; hold on;
% trisurf(tri,x,y,fh(x,y)); 
% view(3);
% title('Exact');
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pu, jac, hess]=JacobiTrig(N, a, b, x, y)

% Get the Jacobi polynomials and their derivatives on a triangle
% 
%  Calling sequence: 
% 
%    Pu=JacobiTrig(N, a, b, x, y)
% 
%    [Pu, jac]=JacobiTrig(N, a, b, x, y)
% 
%    [Pu, jac, hess]=JacobiTrig(N, a, b, x, y)
%   
%  Input :
%
%    N = The index of basis
%
%    a, b - Parameters of Jacobi polynomials
%
%    x, y - A vector of evaluation points on a unit left triangle
% 
%  Output :
% 
%    Pu - The Jacobi polynomials on a unit left triangle
% 
%    jac - The first derivatives of Jacobi polynomials on a unit left triangle
% 
%    hess - The second derivatives of Jacobi polynomials on a unit left triangle
% 

% Coodinate transformation
s=y-x; t=x+y;

% Prepare matrices
n=length(s); 
Pu=zeros(n, N+1); 
if nargout>=2
    jac{1}=Pu; jac{2}=Pu; 
end
if nargout==3
    hess{1,1}=Pu; hess{1,2}=Pu; hess{2,2}=Pu; 
end

% Evaluate
if N==0
    Pu(:,1)=ones(n, 1);
elseif N==1
    P=t*(a-b)/2+s*(a+b+2)/2; 
    Pu(:,1)=ones(n, 1);
    Pu(:,2)=P;
    if nargout>=2
        dPs=(a+b+2)/2; 
        dPt=(a-b)/2; 
        jac{1}(:,2)=(dPt-dPs)*Pu(:,1);
        jac{2}(:,2)=(dPt+dPs)*Pu(:,1);
    end
    if nargout==3
        ddPss=0; 
        ddPst=0; 
        ddPtt=0; 
        hess{1,1}(:,2)=(ddPss-2*ddPst+ddPtt)*Pu(:,1);
        hess{1,2}(:,2)=(ddPtt-ddPss)*Pu(:,1);
        hess{2,2}(:,2)=(ddPss+2*ddPst+ddPtt)*Pu(:,1);
    end
else
    P1=1; P2=t*(a-b)/2+s*(a+b+2)/2;
    Pu(:,1)=ones(n, 1); Pu(:,2)=P2;
    if nargout>=2
        dPs1=0; dPs2=(a+b+2)/2; 
        dPt1=0; dPt2=(a-b)/2; 
        jac{1}(:,2)=(dPt2-dPs2)*Pu(:,1);
        jac{2}(:,2)=(dPt2+dPs2)*Pu(:,1);
    end
    if nargout==3
        ddPss1=0; ddPss2=0; 
        ddPst1=0; ddPst2=0; 
        ddPtt1=0; ddPtt2=0; 
        hess{1,1}(:,2)=(ddPss2-2*ddPst2+ddPtt2)*Pu(:,1);
        hess{1,2}(:,2)=(ddPtt2-ddPss2)*Pu(:,1);
        hess{2,2}(:,2)=(ddPss2+2*ddPst2+ddPtt2)*Pu(:,1);
    end 
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b); 
        a2=(2*n+a+b+1)*(a^2-b^2); 
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2); 
        a4=2*(n+a)*(n+b)*(2*n+a+b+2); 
        P=((a2*t+a3*s).*P2-a4*t.^2.*P1)/a1; 
        Pu(:,n+2)=P; 
        if nargout>=2
            dPs=((a2*t+a3*s).*dPs2+a3*P2-a4*t.^2.*dPs1)/a1; 
            dPt=((a2*t+a3*s).*dPt2+a2*P2-a4*t.^2.*dPt1-2*a4*t.*P1)/a1; 
            jac{1}(:,n+2)=dPt-dPs; 
            jac{2}(:,n+2)=dPt+dPs; 
        end
        if nargout==3
            ddPss=((a2*t+a3*s).*ddPss2+2*a3*dPs2-a4*t.^2.*ddPss1)/a1; 
            ddPst=((a2*t+a3*s).*ddPst2+a2*dPs2+a3*dPt2-2*a4*t.*dPs1-a4*t.^2.*ddPst1)/a1; 
            ddPtt=((a2*t+a3*s).*ddPtt2+2*a2*dPt2-2*a4*P1-4*a4*t.*dPt1-a4*t.^2.*ddPtt1)/a1; 
            hess{1,1}(:,n+2)=ddPss-2*ddPst+ddPtt; 
            hess{1,2}(:,n+2)=ddPtt-ddPss; 
            hess{2,2}(:,n+2)=ddPss+2*ddPst+ddPtt; 
        end
        
        P1=P2; P2=P; 
        if nargout>=2
            dPs1=dPs2; dPs2=dPs; 
            dPt1=dPt2; dPt2=dPt; 
        end
        if nargout==3
            ddPss1=ddPss2; ddPss2=ddPss; 
            ddPst1=ddPst2; ddPst2=ddPst; 
            ddPtt1=ddPtt2; ddPtt2=ddPtt; 
        end
    end
end
if nargout==3
    hess{2,1}=hess{1,2};
end


%% demo
% % Node numbers (n) and the number of basis (N>=0)
% % Parameters of Jacobi polynomials (a, b)
% n=17; N=6; a=2; b=2;
% 
% % The number of points covered by RBF
% nup=3*n; % this parameter can be optimized ...
% 
% % Get nodes
% [x, y]=TrigNodeVect(n);
% tri=tridelaunay(n);
% pnts=[x(:)'; y(:)']; 
% 
% % Get the basis on triangle
% [Pu, jac, hess]=JacobiTrig(N, a, b, x, y);
% dPux=jac{1}; dPuy=jac{2};
% ddPuxx=hess{1,1}; ddPuxy=hess{1,2}; 
% ddPuyy=hess{2,2};
% 
% % Interpolate the first derivetives by RBF
% [G, Dx, Dy]=rbfbasis2d(pnts, pnts, 3*n);
% dPuxi=Dx*Pu;
% dPuyi=Dy*Pu;
% ddPuxxi=Dx*dPux;
% ddPuyxi=Dx*dPuy;
% ddPuxyi=ddPuyxi;
% ddPuyyi=Dy*dPuy;
% 
% k=5;
% figure; hold on;
% trisurf(tri,x,y,dPuy(:,k)); 
% shading interp;
% view(3);
% title('Direct');
% 
% figure; hold on;
% trisurf(tri,x,y,dPuyi(:,k)); 
% shading interp;
% view(3);
% title('RBF');
% 
% figure; hold on;
% trisurf(tri,x,y,dPuy(:,k)-dPuyi(:,k)); 
% shading interp;
% view(3);
% title('Error');
% 
% figure; hold on;
% s=y-x; t=x+y;
% triplot(tri, s,t);
% triplot(tri, s,t,'ro');
% for i=1:length(s)
%     text(s(i), t(i), num2str(i));
% end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dHJ, dJ]=HierarJC0(N, x, der, a, b, m)

% Get the C0 Jacobi hierarchical basis and their derivatives
% 
%  Calling sequence: 
% 
%   [dHJ, dJ]=HierarJC0(N, x, der, a, b, m)
%   
%  Input :
%
%    N = the number of basis
%    x - a vector of evaluation points
%    der - the order of derivatives
%    a, b - parameters of Jacobi polynomials
%    m - the power of the polynomials of coefficients
% 
%  Output :
% 
%    dHL - the C0 Jacobi hierarchical basis and their derivatives
%    dJ - the Jacobi polynomials and their derivatives
% 

if nargin==3
    a=2; b=2; m=1;
end

dHJ=cell(N, der+1);
dJ=JacobiDerCellA(N-1, a, b, x, der);
if m<=1
    p0=(x.^2-1); p1=2*x; p2=2;
else
    p0=(x.^2-1).^m; 
    p1=(m*(x.^2-1).^(m-1)).*(2*x); 
    p2=(m*(m-1)*(x.^2-1).^(m-2)).*(2*x).^2+2*m*(x.^2-1).^(m-1); 
end
if der>=0
    for k=1:N
        dHJ{k, 1}=p0.*dJ{k, 1};
    end
end
if der>=1
    for k=1:N
        dHJ{k, 2}=p1.*dJ{k, 1}+p0.*dJ{k, 2};
    end
end
tp=1;
for t=2:der
    for k=1:N
        dHJ{k, t+1}=tp*p2.*dJ{k, t-1}+t*p1.*dJ{k, t}+p0.*dJ{k, t+1};
    end
    tp=tp+t;
end
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
