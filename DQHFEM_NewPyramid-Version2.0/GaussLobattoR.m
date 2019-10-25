% Solve roots and weights of Gauss-Lobatto quadrature

function [x,C]=GaussLobattoR(n,a,b,R)

if nargin==3
    R=5;
end
rdLp=dLegendreRt(n-1,R);
x=[-1;rdLp;1];

vLp=LegendreRecursion(n-1,rdLp);
C=(2/(n*(n-1)))*[1;1./(vLp.^2);1];

x=(b-a)*x/2+(b+a)/2;
C=(b-a)*C/2;

end

% Get the roots of first order derivative of Legendre polynomials
%
% rxd :    the roots of first order derivative of Legendre polynomials
% m :  the order of Legendre polynomials

function rxd=dLegendreRt(m,R)

    if nargin==1
        R=5;
    end

    n=7*m;
    x=LobattoChebyshev(-1,1,n);
    [~,dy,~]=LegendreRecDer(m,x);

    % Get the span of roots
    t=1;
    rdy=zeros(m-1,2); xdr=rdy;
    nid=rdy(:,1);
    for i=1:fix(n/2)
        if dy(i)*dy(i+1)<0
            rdy(t,1)=dy(i);
            rdy(t,2)=dy(i+1);
            xdr(t,1)=x(i);
            xdr(t,2)=x(i+1);
            nid(t)=i;
            t=t+1;
        end
    end

    % Solve the roots of Legendre
    rxd=zeros(m-1,1); 
    for i=1:fix(m/2)
        if abs(rdy(i,1))<abs(rdy(i,2))
            rxd(i)=xdr(i,1);
        else
            rxd(i)=xdr(i,2);
        end
        for j=1:R
            [~,dyj,ddyj]=LegendreRecDer(m,rxd(i));
            rxd(i)=rxd(i)-dyj/ddyj;
        end
    end

    % Expand
    if mod(m,2)==0
        p=fix(m/2);
        for i=1:p-1
            rxd(p+i)=-rxd(p-i);
        end
    elseif mod(m,2)==1
        p=fix(m/2);
        for i=1:fix(m/2)
            rxd(p+i)=-rxd(p-i+1);
        end
    end
end

% Legendre by the recursion formula

function y=LegendreRecursion(n,x)

    if n==0
        y=1;
    elseif n==1
        y=x;
    elseif n>1    
        y1=1; y2=x;
        for k=1:n-1
            y=((2*k+1)/(k+1))*x.*y2-(k/(k+1))*y1;
            y1=y2;
            y2=y;
        end
    end
end

% Get derivatives of Legendre polynomial by the recursion formula

function [y,dy,ddy]=LegendreRecDer(n,x)

    if n==0
        y=1; dy=0; ddy=0;
    elseif n==1
        y=x; dy=1; ddy=0;
    elseif n>1    
        y1=1; y2=x;
        dy1=1;
        ddy1=0;
        for k=1:n-1
            y=((2*k+1)/(k+1))*x.*y2-(k/(k+1))*y1;
            y1=y2;
            y2=y;

            dy=x.*dy1+(k+1)*y1;
            dy0=dy1;
            dy1=dy;

            ddy=x.*ddy1+(k+2)*dy0;
            ddy1=ddy;
        end
    end
end

% Lobatto Chebyshev nodes
function x=LobattoChebyshev(a,b,N)

    x=zeros(N,1); % 求节点坐标
    for j=1:N
        nd=(j-1)/(N-1); h=0.5*(1-cos(nd*pi));
        x(j)=a+h*(b-a);
    end
end


