function [Nbase,Nnode]=CountBasenNum(number,type)
% 计算number阶次下，单元基函数的总个数
% type：单元类型：1-六面体；2-三棱柱；3-四面体；4-金字塔
switch type
    case 1
        ne=sum(number.edge);
        nf=0;
        for i=1:length(number.face)
            nf=nf+number.face{i}(1)*number.face{i}(2);
        end
        nh=number.H(1)*number.H(2)*number.H(3);
        Nnode=8+ne+nf;
        Nbase=Nnode+nh;
    case 2
        ne=sum(number.edge);
        nf=0;
        for i=2:length(number.face)-1
            nf=nf+number.face{i}(1)*number.face{i}(2);
        end
        nf=nf+number.face{1}*(number.face{1}+1)/2+number.face{5}*(number.face{5}+1)/2;
        nh=number.H(1)*(number.H(1)+1)/2*number.H(2);
        Nnode=6+ne+nf;
        Nbase=Nnode+nh;
    case 3
        ne=sum(number.edge);
        nf=0;
        for i=1:length(number.face)
            nf=nf+number.face{i}*(number.face{i}+1)/2;
        end
        nh=0;
        for i=1:number.H;
            nh=nh+i*(i+1)/2;
        end
        Nnode=4+ne+nf;
        Nbase=Nnode+nh;
    case 4
        ne=sum(number.edge);
        nf=0;
        for i=2:length(number.face)
            nf=nf+number.face{i}*(number.face{i}+1)/2;
        end
        nf=nf+number.face{1}(1)*number.face{1}(2);
        nh=0;
        for i=1:number.H;
            nh=nh+i*(i+1)/2;
        end
        Nnode=5+ne+nf;
        Nbase=Nnode+nh;
end
end
