function NewElementData=VectorRecovery(u,bc,ElementData)
% ���յõ�������װ����ʩ���˱߽�������ȱʡ��ĳЩ�ڵ����λ�ƣ���λ��ʸ����
% ��Ҫ�ָ��õ�ÿ����Ԫ��λ��ʸ�������ܼ���ÿ����Ԫ��λ�Ƴ�
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



