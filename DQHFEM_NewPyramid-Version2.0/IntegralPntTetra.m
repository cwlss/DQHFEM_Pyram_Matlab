function IP=IntegralPntTetra(nq)
% 四面体域积分点
[x,Cx]=GaussLobattoR(nq(1),0,1);
[y,Cy]=GaussLobattoR(nq(2),0,1);
[z,Cz]=GaussLobattoR(nq(3),0,1);
R=zeros(nq(1)*nq(2)*nq(3),1);S=R;T=R;C=R;
ijk=1;
for k=1:nq(3)
    for j=1:nq(2)
        for i=1:nq(1)
            R(ijk)=(1-y(j))*(1-z(k))*x(i);
            S(ijk)=(1-z(k))*y(j);
            T(ijk)=z(k);
            C(ijk)=(1-y(j))*(1-z(k))*(1-z(k))*Cx(i)*Cy(j)*Cz(k);
            ijk=ijk+1;
        end
    end
end
IP.Ri=R;
IP.Si=S;
IP.Ti=T;
IP.Ci=C;
end