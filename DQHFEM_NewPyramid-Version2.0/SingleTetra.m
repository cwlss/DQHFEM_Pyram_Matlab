function ElementData=SingleTetra(a,b,h,n)
% 模型录入程序
% element 1 Ttra
ele1=struct;
ele1.Type=3;
ele1.MapRight={[0;a;0;0],[0;0;b;0],[0;0;0;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele1.Number=number;
ElementData={ele1};
end