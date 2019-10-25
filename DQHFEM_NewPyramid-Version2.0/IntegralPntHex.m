function IP=IntegralPntHex(nq)
% 计算六面体域积分点
[x,Cx]=GaussLobattoR(nq(1),0,1);
[y,Cy]=GaussLobattoR(nq(2),0,1);
[z,Cz]=GaussLobattoR(nq(3),0,1);
R=zeros(nq(1)*nq(2)*nq(3),1);S=R;T=R;C=R;
ijk=1;
for k=1:nq(3)
    for j=1:nq(2)
        for i=1:nq(1)
            R(ijk)=x(i);
            S(ijk)=y(j);
            T(ijk)=z(k);
            C(ijk)=Cx(i)*Cy(j)*Cz(k);
            ijk=ijk+1;
        end
    end
end
IP.Ri=R;
IP.Si=S;
IP.Ti=T;
IP.Ci=C;
end