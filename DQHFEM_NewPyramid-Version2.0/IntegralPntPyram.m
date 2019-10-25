function IP=IntegralPntPyram(nnn)
% 金字塔域积分点
[r,Cr]=GaussLobattoR(nnn(1),-1,1);
[s,Cs]=GaussLobattoR(nnn(2),-1,1);
[t,Ct]=GaussLobattoR(nnn(3),0,1);
N=nnn(1)*nnn(2)*nnn(3);
R=zeros(N,1);S=R;T=R;C=R;
ijk=1;
for i=1:nnn(1)
    for j=1:nnn(2)
        for k=1:nnn(3)
            R(ijk)=r(i)*(1-t(k));
            S(ijk)=s(j)*(1-t(k));
            T(ijk)=t(k);
            C(ijk)=Cr(i)*Cs(j)*Ct(k)*(1-t(k))^2;
            ijk=ijk+1;
        end
    end
end
IP.Ri=R;IP.Si=S;IP.Ti=T;IP.Ci=C;
end