function ElementData=SingleHex(a,b,h,n)
% 模型录入程序
type=1;

ed=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
fa={[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1]};
hh=[n-1,n-1,n-1];
number=struct;
number.edge=ed;number.face=fa;number.H=hh;

ele{1}=[0;a;a;0;0;a;a;0];
ele{2}=[0;0;b;b;0;0;b;b];
ele{3}=[0;0;0;0;h;h;h;h];

Element=struct;
Element.Type=type;
Element.MapRight=ele;
Element.Number=number;
ElementData={Element};

