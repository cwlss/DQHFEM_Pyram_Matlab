function [GK,GM,GP]=Assemble(ElementData)
%% ��ȡȫ�־��� K M����Ԫ�նȾ�����װ��
%GK���ܸ�
%GM����������
%GP�������غ�����
%��ȱʡ����[GK,GM]=Assemble(ElementData)

gn=ElementData{1}.NG;%�������ɶȵĳ���NG�������ڵ�һ����Ԫ�ṹ��ġ�NG���ֶ��У����൥Ԫ�ṹ���ǲ�����NG���ֶ�

Nele=length(ElementData);
GNuvw=cell(Nele,1);
for ni=1:Nele
    GNuvw{ni}=[ElementData{ni}.GlobIndex,ElementData{ni}.GlobIndex+gn,ElementData{ni}.GlobIndex+2*gn]';
end
GK=zeros(3*gn);GP=zeros(3*gn,1);
GM=GK;
for n=1:Nele
    Nb=ElementData{n}.Nbase;
    for i=1:3*Nb
        for j=1:3*Nb
            GK(GNuvw{n}(i),GNuvw{n}(j))=GK(GNuvw{n}(i),GNuvw{n}(j))+ElementData{n}.K(i,j);
            GM(GNuvw{n}(i),GNuvw{n}(j))=GM(GNuvw{n}(i),GNuvw{n}(j))+ElementData{n}.M(i,j);
        end
    end
    if nargout>2
        for k=1:3*Nb
            GP(GNuvw{n}(k))=GP(GNuvw{n}(k))+ElementData{n}.P(k);
        end
    end
end
end


