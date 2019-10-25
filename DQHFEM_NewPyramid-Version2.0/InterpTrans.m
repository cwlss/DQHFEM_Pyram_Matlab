function [GG,GR,GS,GT]=InterpTrans(G,Gr,Gs,Gt,T)
% 使用T进行对G,Gr,Gs,Gt转换
GG=G*T;
GR=Gr*T;
GS=Gs*T;
GT=Gt*T;
end