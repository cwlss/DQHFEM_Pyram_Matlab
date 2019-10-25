function [X,Y,Z,dX,dY,dZ,dJ]=MapPrism(r,s,t,x,y,z)
nr=length(r);
g1=(1-r-s).*(1-t);
g2=r.*(1-t);
g3=s.*(1-t);
g4=(1-r-s).*t;
g5=r.*t;
g6=s.*t;
g=[g1,g2,g3,g4,g5,g6];
X=g*x;Y=g*y;Z=g*z;
if nargout>3
    dg1dr=-(1-t);   dg1ds=-(1-t);    dg1dt=-(1-r-s);
    dg2dr=(1-t);    dg2ds=0*s;       dg2dt=-r;
    dg3dr=0*r;      dg3ds=(1-t);     dg3dt=-s;
    dg4dr=-t;       dg4ds=-t;        dg4dt=(1-r-s);
    dg5dr=t;        dg5ds=0*s;       dg5dt=r;
    dg6dr=0*r;      dg6ds=t;         dg6dt=s;
    
    dgdr=[dg1dr,dg2dr,dg3dr,dg4dr,dg5dr,dg6dr];
    dgds=[dg1ds,dg2ds,dg3ds,dg4ds,dg5ds,dg6ds];
    dgdt=[dg1dt,dg2dt,dg3dt,dg4dt,dg5dt,dg6dt];
    
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