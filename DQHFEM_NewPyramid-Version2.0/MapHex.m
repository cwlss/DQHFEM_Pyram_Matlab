function [X,Y,Z,dX,dY,dZ,dJ]=MapHex(r,s,t,x,y,z)
%直边映射，可被MapCurHex替代
nr=length(r);
g1=(1-r).*(1-s).*(1-t);
g2=r.*(1-s).*(1-t);
g3=r.*s.*(1-t);
g4=(1-r).*s.*(1-t);
g5=(1-r).*(1-s).*t;
g6=r.*(1-s).*t;
g7=r.*s.*t;
g8=(1-r).*s.*t;
g=[g1,g2,g3,g4,g5,g6,g7,g8];
X=g*x;Y=g*y;Z=g*z;
if nargout>3
    dg1dr=-(1-s).*(1-t);   dg1ds=-(1-r).*(1-t);   dg1dt=-(1-r).*(1-s);
    dg2dr=(1-s).*(1-t);    dg2ds=-r.*(1-t);       dg2dt=-r.*(1-s);
    dg3dr=s.*(1-t);        dg3ds=r.*(1-t);        dg3dt=-r.*s;
    dg4dr=-s.*(1-t);       dg4ds=(1-r).*(1-t);    dg4dt=-(1-r).*s;
    dg5dr=-(1-s).*t;       dg5ds=-(1-r).*t;       dg5dt=(1-r).*(1-s);
    dg6dr=(1-s).*t;        dg6ds=-r.*t;           dg6dt=r.*(1-s);
    dg7dr=s.*t;            dg7ds=r.*t;            dg7dt=r.*s;
    dg8dr=-s.*t;           dg8ds=(1-r).*t;        dg8dt=(1-r).*s;
    dgdr=[dg1dr,dg2dr,dg3dr,dg4dr,dg5dr,dg6dr,dg7dr,dg8dr];
    dgds=[dg1ds,dg2ds,dg3ds,dg4ds,dg5ds,dg6ds,dg7ds,dg8ds];
    dgdt=[dg1dt,dg2dt,dg3dt,dg4dt,dg5dt,dg6dt,dg7dt,dg8dt];
    
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