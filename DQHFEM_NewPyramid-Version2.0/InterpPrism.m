function T=InterpPrism(number,opt)
% �������������׻���������ֵ��������ת������T���������ֵת��
% opt��1 or ȱʡ-Fekete���ֵ��  2-���Ƚڵ��ֵ
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
    