function [X,Y,Z,dX,dY,dZ,dJ]=MapTetra(r,s,t,x,y,z)
nr=length(r);
g1=1-r-s-t;
g2=r;
g3=s;
g4=t;

g=[g1,g2,g3,g4];
X=g*x;Y=g*y;Z=g*z;
if nargout>3
    I=ones(size(r));
    dg1dr=-I;   dg1ds=-I;    dg1dt=-I;
    dg2dr=I;    dg2ds=0*s;       dg2dt=0*t;
    dg3dr=0*r;      dg3ds=I;     dg3dt=0*t;
    dg4dr=0*r;       dg4ds=0*s;        dg4dt=I;

    dgdr=[dg1dr,dg2dr,dg3dr,dg4dr];
    dgds=[dg1ds,dg2ds,dg3ds,dg4ds];
    dgdt=[dg1dt,dg2dt,dg3dt,dg4dt];
    
    dxdr=dgdr*x;dxds=dgds*x;dxdt=dgdt*x;
    dX={dxdr;dxds;dxdt};
    
    dydr=dgdr*y;dyds=dgds*y;dydt=dgdt*y;
    dY={dydr;dyds;dydt};
    
    dzdr=dgdr*z;dzds=dgds*z;dzdt=dgdt*z;
    dZ={dzdr;dzds;dzdt};
    dJ=zeros(nr,1);
    for i=1:nr
        J=[dxdr(i),dydr(i),dzdr(i)
            dxds(i),dyds(i),dzds(i)
            dxdt(i),dydt(i),dzdt(i)];
        dJ(i)=det(J);
    end
end
end
%%test

%S=s;R=r;T=t;
% Jx=1/dJ*[dyds*dzdt-dydt*dzds,dydt*dzdr-dydr*dzdt,dydr*dzds-dyds*dzdr];
% Jy=1/dJ*[dxdt*dzds-dxds*dzdt,dxdr*dzdt-dxdt*dzdr,dxds*dzdr-dxdr*dzds];
% Jz=1/dJ*[dxds*dydt-dxdt*dyds,dxdt*dydr-dxdr*dydt,dxdr*dyds-dxds*dydr];
% if nargout>4
%     X=x*g;
%     Y=y*g;
%     Z=z*g;
% end
% end