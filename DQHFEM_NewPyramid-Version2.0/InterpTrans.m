function [GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T)
% ʹ��T���ж�G,Gr,Gs,Gtת��
GG=G*T;
GR=Gr*T;
GS=Gs*T;
GT=Gt*T;
end