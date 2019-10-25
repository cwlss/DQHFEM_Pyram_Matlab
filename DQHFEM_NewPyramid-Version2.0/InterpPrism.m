function T=InterpPrism(number,opt)
% 计算三棱柱阶谱基函数到插值基函数的转换矩阵T：六面体插值转换
% opt：1 or 缺省-Fekete点插值；  2-均匀节点插值
if nargin==1
[Rb,Sb,Tb]=FeketePrism(number);
else if opt==1
        [Rb,Sb,Tb]=FeketePrism(number);
    else if opt==2
        [Rb,Sb,Tb]=UniformPrism(number);
        else
            error('error')
        end
    end
end
g=BasePrism(Rb,Sb,Tb,number);
[N1,N2]=size(g);
G1=g(:,1:N1);
G1_=G1\eye(N1);
I=eye(N2-N1);
O1=zeros(N1,N2-N1);
O2=zeros(N2-N1,N1);
T=[G1_,O1
    O2,I];
end
    