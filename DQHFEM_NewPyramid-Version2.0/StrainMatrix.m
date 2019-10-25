function B=StrainMatrix(Gx,Gy,Gz)
O=zeros(size(Gx));
B=[Gx,O,O
    O,Gy,O
    O,O,Gz
    Gy,Gx,O
    O,Gz,Gy
    Gz,O,Gx];
end
    