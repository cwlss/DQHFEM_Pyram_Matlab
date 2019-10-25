function ElementData=DoubleTetraPyramid2(a,b,h,n)
%% 模型录入程序
%%
a=2*a;
% element 1 Ttra
ele1=struct;
ele1.Type=3;
ele1.MapRight={[0;a/(2^0.5);0;0],[0;0;a/(2^0.5);0],[0;0;0;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele1.Number=number;


% element 2 Ttra
ele2=struct;
ele2.Type=3;
ele2.MapRight={[0;0;a/(2^0.5);0],[0;-a/(2^0.5);0;0],[0;0;0;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele2.Number=number;
ElementData={ele1,ele2};
end