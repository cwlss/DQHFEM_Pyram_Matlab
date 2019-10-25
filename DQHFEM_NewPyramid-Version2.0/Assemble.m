function [GK,GM,GP]=Assemble(ElementData)
%% 获取全局矩阵 K M（单元刚度矩阵组装）
%GK：总刚
%GM：质量矩阵
%GP：总体载荷列阵
%可缺省调用[GK,GM]=Assemble(ElementData)

gn=ElementData{1}.NG;%总体自由度的长度NG，储存在第一个单元结构体的‘NG’字段中，其余单元结构体是不含‘NG’字段

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


