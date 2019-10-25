function NewElementData=VectorRecovery(u,bc,ElementData)
% 最终得到的是组装后，且施加了边界条件（缺省了某些节点的零位移）的位移矢量，
% 需要恢复得到每个单元的位移矢量，才能计算每个单元的位移场
Nele=length(ElementData);
N=length(bc);
L=N/3;
bn=find(bc);
U=zeros(N,1);
for i=1:length(bn)
    U(bn(i))=u(i);
end
for i=1:Nele
    ElementData{i}.eleU=U(ElementData{i}.GlobIndex);
    ElementData{i}.eleV=U(ElementData{i}.GlobIndex+L);
    ElementData{i}.eleW=U(ElementData{i}.GlobIndex+2*L);

end
NewElementData=ElementData;
end



