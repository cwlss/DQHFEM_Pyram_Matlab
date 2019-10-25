function ElementData=SinglePyramid(a,b,h,n)
% 模型录入程序
type=4;

ed=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
fa={[n-1,n-1],n-2,n-2,n-2,n-2};
hh=n-2;
number=struct;number.edge=ed;number.face=fa;number.H=hh;
%%
MapRight{1}=[-a; a; a;-a; 0];
MapRight{2}=[-b;-b; b; b; 0];
MapRight{3}=[0 ; 0; 0; 0; h];

Element=struct;
Element.Type=type;
Element.MapRight=MapRight;
Element.Number=number;
ElementData={Element};