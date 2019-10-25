
function DrawBasePyram(m,BaseType,opt1,opt2,NodeSwitch)
% Input:
% BaseType:可视化基函数的类型（1——顶点基；2——边基；3——面基；4——体基）
% opt1:可视化基函数的导数类型（0——函数值G；1——Gr；2——Gs；3——Gt）
% opt2:0 or缺省——非插值形函数，1——Fekete插值形函数
% m:基函数阶次选择，m建议别选太大，否则基函数数量太多，运行时长太久
% 体函数由于数量较多，只显示前10个基

number.edge=[m-1,m-1,m-1,m-1,m-1,m-1,m-1,m-1];
number.face={[m-1,m-1],m-2,m-2,m-2,m-2};
number.H=m-2;
n=50;
%% 金字塔网格点
% [r,s,t,tri]=StdPyramPoint(n);
[R, S, T, ~]=meshgridpyramc(n, n, n);
r=R(:);s=S(:);t=T(:);
[G,Gr,Gs,Gt]=BasePyram(r,s,t,number);
if nargin>3 && opt2==1
    Tr=InterpPyram(number);
    [G,Gr,Gs,Gt]=InterpTrans(G,Gr,Gs,Gt,Tr);
end

switch BaseType
    case 1
        for i=1:5
            set(0,'defaultfigurecolor','w');
            figure; hold on;
            title(['Vertex',num2str(i)]);
            
            switch opt1
                case 0
                    c=G(:,i);
                case 1
                    c=Gr(:,i);
                case 2
                    c=Gs(:,i);
                case 3
                    c=Gt(:,i);
            end
            C=reshape(c, n, n, n);
            hexasurfg(R, S, T, C, 2);
            
            colorbar;
            view(3); alpha(1);
            axis equal;
            axis off
            shading interp
            if nargin>4&&NodeSwitch==1
                [Rb,Sb,Tb]=FeketePyram(number);
                scatter3(Rb,Sb,Tb,70,'filled','r');
            end
            hold on
            
            
        end
    case 2
        for i=6:5+8*(m-1)
            set(0,'defaultfigurecolor','w');
            figure; hold on;
            title(['Edge',num2str(i-5)])
            
            switch opt1
                case 0
                    c=G(:,i);
                case 1
                    c=Gr(:,i);
                case 2
                    c=Gs(:,i);
                case 3
                    c=Gt(:,i);
            end
            C=reshape(c, n, n, n);
            hexasurfg(R, S, T, C, 2);
            
            colorbar;
            view(3);
            alpha(1);
            axis equal;
            axis off
            shading interp
            if nargin>4&&NodeSwitch==1
                [Rb,Sb,Tb]=FeketePyram(number);
                scatter3(Rb,Sb,Tb,70,'filled','r');
            end
            hold on
        end
    case 3
        for i=6+8*(m-1):5+8*(m-1)+(m-1)^2+2*(m-1)*(m-2)
            set(0,'defaultfigurecolor','w');
            figure; hold on;
            title(['Face',num2str(i-6-8*(m-1))])
            switch opt1
                case 0
                    c=G(:,i);
                case 1
                    c=Gr(:,i);
                case 2
                    c=Gs(:,i);
                case 3
                    c=Gt(:,i);
            end
            C=reshape(c, n, n, n);
            hexasurfg(R, S, T, C, 2);
            colorbar;
            view(3); alpha(1);
            axis equal;
            axis off
            shading interp
            if nargin>4&&NodeSwitch==1
                [Rb,Sb,Tb]=FeketePyram(number);
                scatter3(Rb,Sb,Tb,70,'filled','r');
            end
            hold on
        end
    case 4
        for i=6+8*(m-1)+(m-1)^2+2*(m-1)*(m-2):5+8*(m-1)+(m-1)^2+2*(m-1)*(m-2)+10
            set(0,'defaultfigurecolor','w');
            figure; hold on;
            title(['Body',num2str(i-5-8*(m-1)-(m-1)^2-2*(m-1)*(m-2))])
            switch opt1
                case 0
                    c=G(:,i);
                case 1
                    c=Gr(:,i);
                case 2
                    c=Gs(:,i);
                case 3
                    c=Gt(:,i);
            end
            C=reshape(c, n, n, n);
            hexasurfg(R, S, T, C, 2);
            hold on
            colorbar;
            view(3); alpha(0.1);
            axis equal;
            axis off
            shading interp
            if nargin>4&&NodeSwitch==1
                [Rb,Sb,Tb]=FeketePyram(number);
                scatter3(Rb,Sb,Tb,70,'filled','r');
            end
            hold on
        end
end

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

function hexasurfg(x, y, z, C, k)
%
% hexasurfg : Plot a hexahedron region by surfaces
%
% Calling Sequences:
%
%     hexasurfg(x, y, z)
%
%     hexasurfg(x, y, z, C)
%
%     hexasurfg(x, y, z, C, k)
%
% INPUTS:
%
%     x, y, z   -    Matrices of parametric coordinates of a unit
%               hexahedron.  See also  meshgridtetrac.
%
%     C - A matrix of color
%
%     k  -  Step leangth of plot. For N grid, the surfaces 1:k:N
%            of each direction will be plotted.
%

if nargin==3
    C=z; k=1;
elseif nargin==4
    k=1;
end

hold_flag = ishold;
hold on;

M=size(y, 1); N=size(y, 2); K=size(y, 3);
for i=1:k:M
    surf(squeeze(x(i,:,:)), squeeze(y(i,:,:)), squeeze(z(i,:,:)), squeeze(C(i,:,:)));
end
for i=1:k:N
    surf(squeeze(x(:,i,:)), squeeze(y(:,i,:)), squeeze(z(:,i,:)), squeeze(C(:,i,:)));
end
for i=1:k:K
    surf(squeeze(x(:,:,i)), squeeze(y(:,:,i)), squeeze(z(:,:,i)), squeeze(C(:,:,i)));
end
shading interp;
colorbar;
view(3); alpha(0.1);

if (~hold_flag)
    hold off
end
end

