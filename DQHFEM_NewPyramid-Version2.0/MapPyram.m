function [X,Y,Z,dX,dY,dZ,dJ]=MapPyram(r,s,t,x,y,z)
%
nr=length(r);
g1=zeros(size(r));g2=g1;g3=g1;g4=g1;g5=g1;
dg1dr=g1;dg2dr=g1;dg3dr=g1;dg4dr=g1;dg5dr=g1;
dg1ds=g1;dg2ds=g1;dg3ds=g1;dg4ds=g1;dg5ds=g1;
dg1dt=g1;dg2dt=g1;dg3dt=g1;dg4dt=g1;dg5dt=g1;
for i=1:nr
    if t(i)==1
        g1(i)=0;
        g2(i)=0;
        g3(i)=0;
        g4(i)=0;
        g5(i)=1;
    else
        g1(i)=0.25*(1-r(i)-s(i)-t(i)+r(i)*s(i)/(1-t(i)));
        g2(i)=0.25*(1+r(i)-s(i)-t(i)-r(i)*s(i)/(1-t(i)));
        g3(i)=0.25*(1+r(i)+s(i)-t(i)+r(i)*s(i)/(1-t(i)));
        g4(i)=0.25*(1-r(i)+s(i)-t(i)-r(i)*s(i)/(1-t(i)));
        g5(i)=t(i);
    end
    %%
    if nargout>3
        if t(i)==1
            dg1dr(i)=0.25*(-1);   dg1ds(i)=0.25*(-1);   dg1dt(i)=0.25.*(-1);
            dg2dr(i)=0.25.*(1);   dg2ds(i)=0.25.*(-1);  dg2dt(i)=0.25.*(-1);
            dg3dr(i)=0.25.*(1);   dg3ds(i)=0.25.*(1);   dg3dt(i)=0.25.*(-1);
            dg4dr(i)=0.25.*(-1);  dg4ds(i)=0.25.*(1);   dg4dt(i)=0.25.*(-1);
            dg5dr(i)=0;           dg5ds(i)=0;           dg5dt(i)=1;
        else
            dg1dr(i)=0.25*(-1+s(i)/(1-t(i)));   dg1ds(i)=0.25*(-1+r(i)/(1-t(i)));   dg1dt(i)=0.25*(-1+r(i)*s(i)/(1-t(i)).^2);
            dg2dr(i)=0.25*(1-s(i)/(1-t(i)));    dg2ds(i)=0.25*(-1-r(i)/(1-t(i)));   dg2dt(i)=0.25*(-1-r(i)*s(i)/(1-t(i)).^2);
            dg3dr(i)=0.25*(1+s(i)/(1-t(i)));    dg3ds(i)=0.25*(1+r(i)/(1-t(i)));    dg3dt(i)=0.25*(-1+r(i)*s(i)/(1-t(i)).^2);
            dg4dr(i)=0.25*(-1-s(i)/(1-t(i)));   dg4ds(i)=0.25*(1-r(i)/(1-t(i)));    dg4dt(i)=0.25*(-1-r(i)*s(i)/(1-t(i)).^2);
            dg5dr(i)=0;                    dg5ds(i)=0;                    dg5dt(i)=1;
        end
    end
    %%
end
g=[g1,g2,g3,g4,g5];
X=g*x;Y=g*y;Z=g*z;
%%
if nargout>3
    dgdr=[dg1dr,dg2dr,dg3dr,dg4dr,dg5dr];
    dgds=[dg1ds,dg2ds,dg3ds,dg4ds,dg5ds];
    dgdt=[dg1dt,dg2dt,dg3dt,dg4dt,dg5dt];
    
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
%%
end

