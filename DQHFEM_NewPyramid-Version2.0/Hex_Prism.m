function ElementData=Hex_Prism(a,b,h,n)
%% 模型录入程序
%%
type=1;
ed=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
fa={[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1]};
hh=[n-1,n-1,n-1];
number=struct;
number.edge=ed;number.face=fa;number.H=hh;

ele{1}=[0;a/2;a/2;0;0;a/2;a/2;0];
ele{2}=[0;0;b;b;0;0;b;b];
ele{3}=[0;0;0;0;h;h;h;h];

Element=struct;
Element.Type=type;
Element.MapRight=ele;
Element.Number=number;

% element2 Prism
ele2=struct;
ele2.Type=2;
ele2.MapRight={[a/2;a;a/2;a/2;a;a/2],[0;0;b;0;0;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele2.Number=number;

% element3 Prism
ele3=struct;
ele3.Type=2;
ele3.MapRight={[a;a;a/2;a;a;a/2],[0;b;b;0;b;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele3.Number=number;

ElementData={Element;ele2;ele3};
%  ElementData={ele2;ele3};