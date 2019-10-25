function K=StiffnessMatrix(B,D,J,C)
JC=J.*C;
JC=diag(JC);
D=kron(D,JC);
K=B'*D*B;
end