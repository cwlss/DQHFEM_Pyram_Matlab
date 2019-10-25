function M=MassMatrix(G,J,C,rhom)
JC=J.*C;
JC=diag(JC);
m=G'*JC*G;
O=zeros(size(m));
M=rhom*[m,O,O
        O,m,O
        O,O,m];
end

