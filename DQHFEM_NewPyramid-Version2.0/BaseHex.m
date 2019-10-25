function [G,Gr,Gs,Gt]=BaseHex(r,s,t,number )
%% 求取六面体基的函数值和三个偏导数值
% number为结构体：储存六面体单元边、面、体上的阶次
% 可缺省调用[G]=BaseHex(r,s,t,number )
%   此处显示详细说明
nedge=number.edge;
nface=number.face;
nH=number.H;
if nargout==1
    S=BaseVertex(r,s,t);
    G=S;
    for i=1:12
        E=BaseEdge(r,s,t,nedge(i),i);
        G=[G,E];
    end
    for j=1:6
        F=BaseFace(r,s,t,nface{j},j);
        G=[G,F];
    end
    H=BaseInterior(r,s,t,nH);
    G=[G,H];
else
    [S,Sr,Ss,St]=BaseVertex(r,s,t);
    G=S;Gr=Sr;Gs=Ss;Gt=St;
    for i=1:12
        [E,Er,Es,Et]=BaseEdge(r,s,t,nedge(i),i);
        G=[G,E];Gr=[Gr,Er];Gs=[Gs,Es];Gt=[Gt,Et];
    end
    for j=1:6
        [F,Fr,Fs,Ft]=BaseFace(r,s,t,nface{j},j);
        G=[G,F];Gr=[Gr,Fr];Gs=[Gs,Fs];Gt=[Gt,Ft];
    end
    [H,Hr,Hs,Ht]=BaseInterior(r,s,t,nH);
    G=[G,H];Gr=[Gr,Hr];Gs=[Gs,Hs];Gt=[Gt,Ht];
end
end


function [H,Hr,Hs,Ht]=BaseInterior(r,s,t,n)
if length(n)~=3
    error('n input is error');
end
nr=length(r);
pe=[1,1,1,1,1,1];
x=2*r-1;y=2*s-1;z=2*t-1;
if nargout==1
    P=BasePolyPart(r,s,t,pe);
    J1=BaseJacobi(n(1)-1, 2, 2, x);
    J2=BaseJacobi(n(2)-1, 2, 2, y);
    J3=BaseJacobi(n(3)-1, 2, 2, z);
    H=zeros(nr,n(1)*n(2)*n(3));
    ijk=1;
    for i=1:n(1)
        for j=1:n(2)
            for k=1:n(3)
                H(:,ijk)=P.*J1(:,i).*J2(:,j).*J3(:,k);
                ijk=ijk+1;
            end
        end
    end
else
    [P,Pr,Ps,Pt]=BasePolyPart(r,s,t,pe);
    [J1,dJ1]=BaseJacobi(n(1)-1, 2, 2, x);
    [J2,dJ2]=BaseJacobi(n(2)-1, 2, 2, y);
    [J3,dJ3]=BaseJacobi(n(3)-1, 2, 2, z);
    H=zeros(nr,n(1)*n(2)*n(3));Hr=H;Hs=H;Ht=H;
    ijk=1;
    for i=1:n(1)
        for j=1:n(2)
            for k=1:n(3)
                H(:,ijk)=P.*J1(:,i).*J2(:,j).*J3(:,k);
                Hr(:,ijk)=Pr.*J1(:,i).*J2(:,j).*J3(:,k)+P.*dJ1(:,i).*J2(:,j).*J3(:,k)*2;
                Hs(:,ijk)=Ps.*J1(:,i).*J2(:,j).*J3(:,k)+P.*J1(:,i).*dJ2(:,j).*J3(:,k)*2;
                Ht(:,ijk)=Pt.*J1(:,i).*J2(:,j).*J3(:,k)+P.*J1(:,i).*J2(:,j).*dJ3(:,k)*2;
                ijk=ijk+1;
            end
        end
    end
end
end

function [F,Fr,Fs,Ft]=BaseFace(r,s,t,n,Nface)
if length(n)~=2
    error('n input is error');
end
nr=length(r);
switch Nface
    case 1
        pe=[1,1,1,1,0,1];
        x=2*r-1;y=2*s-1;
        xr=2;xs=0;xt=0;
        yr=0;ys=2;yt=0;
    case 2
        pe=[1,1,0,1,1,1];
        x=2*r-1;y=2*t-1;
        xr=2;xs=0;xt=0;
        yr=0;ys=0;yt=2;
    case 3
        pe=[1,0,1,1,1,1];
        x=2*s-1;y=2*t-1;
        xr=0;xs=2;xt=0;
        yr=0;ys=0;yt=2;
    case 4
        pe=[1,1,1,0,1,1];
        x=2*r-1;y=2*t-1;
        xr=2;xs=0;xt=0;
        yr=0;ys=0;yt=2;
    case 5
        pe=[0,1,1,1,1,1];
        x=2*s-1;y=2*t-1;
        xr=0;xs=2;xt=0;
        yr=0;ys=0;yt=2;
    case 6
        pe=[1,1,1,1,1,0];
        x=2*r-1;y=2*s-1;
        xr=2;xs=0;xt=0;
        yr=0;ys=2;yt=0;
end
if nargout==1
    J1=BaseJacobi(n(1)-1, 2, 2, x);
    J2=BaseJacobi(n(2)-1, 2, 2, y);
    P=BasePolyPart(r,s,t,pe);
    F=zeros(nr,n(1)*n(2));
    ij=1;
    for i=1:n(1)
        for j=1:n(2)
            F(:,ij)=P.*J1(:,i).*J2(:,j);
            ij=ij+1;
        end
    end
else
    [J1,dJ1]=BaseJacobi(n(1)-1, 2, 2, x);
    [J2,dJ2]=BaseJacobi(n(2)-1, 2, 2, y);
    [P,Pr,Ps,Pt]=BasePolyPart(r,s,t,pe);
    F=zeros(nr,n(1)*n(2));Fr=F;Fs=F;Ft=F;
    ij=1;
    for i=1:n(1)
        for j=1:n(2)
            F(:,ij)=P.*J1(:,i).*J2(:,j);
            Fr(:,ij)=Pr.*J1(:,i).*J2(:,j)+P.*dJ1(:,i).*J2(:,j)*xr+P.*J1(:,i).*dJ2(:,j)*yr;
            Fs(:,ij)=Ps.*J1(:,i).*J2(:,j)+P.*dJ1(:,i).*J2(:,j)*xs+P.*J1(:,i).*dJ2(:,j)*ys;
            Ft(:,ij)=Pt.*J1(:,i).*J2(:,j)+P.*dJ1(:,i).*J2(:,j)*xt+P.*J1(:,i).*dJ2(:,j)*yt;
            ij=ij+1;
        end
    end
end
end

function [E,Er,Es,Et]=BaseEdge(r,s,t,n,Nedge)
nr=length(r);
switch Nedge
    case 1
        pe=[1,1,0,1,0,1];
        x=2*r-1;
        xr=2;xs=0;xt=0;
    case 2
        pe=[1,0,1,1,0,1];
        x=2*s-1;
        xr=0;xs=2;xt=0;
    case 3
        pe=[1,1,1,0,0,1];
        x=2*r-1;
        xr=2;xs=0;xt=0;
    case 4
        pe=[0,1,1,1,0,1];
        x=2*s-1;
        xr=0;xs=2;xt=0;
    case 5
        pe=[0,1,0,1,1,1];
        x=2*t-1;
        xr=0;xs=0;xt=2;
    case 6
        pe=[1,0,0,1,1,1];
        x=2*t-1;
        xr=0;xs=0;xt=2;
    case 7
        pe=[1,0,1,0,1,1];
        x=2*t-1;
        xr=0;xs=0;xt=2;
    case 8
        pe=[0,1,1,0,1,1];
        x=2*t-1;
        xr=0;xs=0;xt=2;
    case 9
        pe=[1,1,0,1,1,0];
        x=2*r-1;
        xr=2;xs=0;xt=0;
    case 10
        pe=[1,0,1,1,1,0];
        x=2*s-1;
        xr=0;xs=2;xt=0;
    case 11
        pe=[1,1,1,0,1,0];
        x=2*r-1;
        xr=2;xs=0;xt=0;
    case 12
        pe=[0,1,1,1,1,0];
        x=2*s-1;
        xr=0;xs=2;xt=0;
end
if nargout==1
    J=BaseJacobi(n-1, 2, 2, x);
    P=BasePolyPart(r,s,t,pe);
    E=zeros(nr,n);
    for i=1:n
        E(:,i)=P.*J(:,i);
    end
else
    [J,dJ]=BaseJacobi(n-1, 2, 2, x);
    [P,Pr,Ps,Pt]=BasePolyPart(r,s,t,pe);
    E=zeros(nr,n);Er=E;Es=E;Et=E;
    for i=1:n
        E(:,i)=P.*J(:,i);
        Er(:,i)=Pr.*J(:,i)+P.*dJ(:,i)*xr;
        Es(:,i)=Ps.*J(:,i)+P.*dJ(:,i)*xs;
        Et(:,i)=Pt.*J(:,i)+P.*dJ(:,i)*xt;
    end
end
end

function [S,Sr,Ss,St]=BaseVertex(r,s,t)
if nargout==1
    S1=(1-r).*(1-s).*(1-t);
    S2=r.*(1-s).*(1-t);
    S3=r.*s.*(1-t);
    S4=(1-r).*s.*(1-t);
    S5=(1-r).*(1-s).*t;
    S6=r.*(1-s).*t;
    S7=r.*s.*t;
    S8=(1-r).*s.*t;
    S=[S1,S2,S3,S4,S5,S6,S7,S8];
else
    S1=(1-r).*(1-s).*(1-t);
    dS1dr=-(1-s).*(1-t);
    dS1ds=-(1-r).*(1-t);
    dS1dt=-(1-r).*(1-s);
    
    S2=r.*(1-s).*(1-t);
    dS2dr=(1-s).*(1-t);
    dS2ds=-r.*(1-t);
    dS2dt=-r.*(1-s);
    
    S3=r.*s.*(1-t);
    dS3dr=s.*(1-t);
    dS3ds=r.*(1-t);
    dS3dt=-r.*s;
    
    S4=(1-r).*s.*(1-t);
    dS4dr=-s.*(1-t);
    dS4ds=(1-r).*(1-t);
    dS4dt=-(1-r).*s;
    
    S5=(1-r).*(1-s).*t;
    dS5dr=-(1-s).*t;
    dS5ds=-(1-r).*t;
    dS5dt=(1-r).*(1-s);
    
    S6=r.*(1-s).*t;
    dS6dr=(1-s).*t;
    dS6ds=-r.*t;
    dS6dt=r.*(1-s);
    
    S7=r.*s.*t;
    dS7dr=s.*t;
    dS7ds=r.*t;
    dS7dt=r.*s;
    
    S8=(1-r).*s.*t;
    dS8dr=-s.*t;
    dS8ds=(1-r).*t;
    dS8dt=(1-r).*s;
    S=[S1,S2,S3,S4,S5,S6,S7,S8];
    Sr=[dS1dr,dS2dr,dS3dr,dS4dr,dS5dr,dS6dr,dS7dr,dS8dr];
    Ss=[dS1ds,dS2ds,dS3ds,dS4ds,dS5ds,dS6ds,dS7ds,dS8ds];
    St=[dS1dt,dS2dt,dS3dt,dS4dt,dS5dt,dS6dt,dS7dt,dS8dt];
end


end

function [P,Pr,Ps,Pt]=BasePolyPart(r,s,t,pe)
A=pe(1);B=pe(2);C=pe(3);D=pe(4);E=pe(5);F=pe(6);
a=r.^A;
b=(1-r).^B;
c=s.^C;
d=(1-s).^D;
e=t.^E;
f=(1-t).^F;
if nargout==1
    P=a.*b.*c.*d.*e.*f;
else
    
    if A==0
        da=0*t;
    else
        da=A*r.^(A-1);
    end
    if B==0
        db=0*t;
    else
        db=-B*(1-r).^(B-1);
    end
    if C==0
        dc=0*t;
    else
        dc=C*s.^(C-1);
    end
    if D==0
        dd=0*t;
    else
        dd=-D*(1-s).^(D-1);
    end
    if E==0
        de=0*t;
    else
        de=E*t.^(E-1);
    end
    if F==0
        df=0*t;
    else
        df=-F*(1-t).^(F-1);
    end
    P=a.*b.*c.*d.*e.*f;
    Pr=c.*d.*e.*f.*(da.*b+a.*db);
    Ps=a.*b.*e.*f.*(dc.*d+c.*dd);
    Pt=a.*b.*c.*d.*(de.*f+e.*df);
end
end

function [J, dJ, ddJ, dddJ]=BaseJacobi(N, a, b, x)
nx=length(x);
J=zeros(nx,N+1);dJ=J;ddJ=J;dddJ=J;
if N==0
    y=ones(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
    J(:,1)=y;dJ(:,1)=dy;ddJ(:,1)=ddy;dddJ(:,1)=dddy;
elseif N==1
    y=ones(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
    J(:,1)=y;dJ(:,1)=dy;ddJ(:,1)=ddy;dddJ(:,1)=dddy;
    y=0.5*(a-b+(a+b+2)*x); dy=0.5*(a+b+2); ddy=0*y; dddy=0*y;
    J(:,2)=y;dJ(:,2)=dy;ddJ(:,2)=ddy;dddJ(:,2)=dddy;
else
    y=ones(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
    J(:,1)=y;dJ(:,1)=dy;ddJ(:,1)=ddy;dddJ(:,1)=dddy;
    y=0.5*(a-b+(a+b+2)*x); dy=0.5*(a+b+2); ddy=0*y; dddy=0*y;
    J(:,2)=y;dJ(:,2)=dy;ddJ(:,2)=ddy;dddJ(:,2)=dddy;
    
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
        J(:,n+2)=y;
        if nargout>1
            dy=((a2+a3*x).*dy2+a3*y2-a4*dy1)/a1;
            dJ(:,n+2)=dy;
        end
        if nargout>2
            ddy=((a2+a3*x).*ddy2+2*a3*dy2-a4*ddy1)/a1;
            ddJ(:,n+2)=ddy;
        end
        if nargout>3
            dddy=((a2+a3*x).*dddy2+3*a3*ddy2-a4*dddy1)/a1;
            dddJ(:,n+2)=dddy;
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

