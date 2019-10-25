function ElementData=ThreePyramidRevised(a,b,h,n)
%% 模型录入程序:优化了厚度方向的节点数目

%% element1 Pyramid
ele1=struct;
ele1.Type=4;
ele1.MapRight={[0;a;a;0;0],[0;0;b;b;0],[0;0;0;0;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,1,n-1,n-1,n-1];
number.face={[n-1,n-1],1,n-2,n-2,1};
number.H=n-2;
ele1.Number=number;

%% element2 Pyramid
ele2=struct;
ele2.Type=4;
ele2.MapRight={[a;a;a;a;0],[0;0;b;b;0],[0;h;h;0;h]};
number=struct;
number.edge=[1,n-1,1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],1,n-2,1,n-2};
number.H=n-2;
ele2.Number=number;
%% element3 Pyramid
ele3=struct;
ele3.Type=4;
ele3.MapRight={[0;a;a;0;0],[b;b;b;b;0],[0;0;h;h;h]};
number=struct;
number.edge=[n-1,1,n-1,1,n-1,n-1,n-1,n-1];
number.face={[n-1,1],n-2,1,n-2,1};
number.H=n-2;
ele3.Number=number;

ElementData={ele1;ele2;ele3};
% ElementData={ele1;ele2;ele3;ele4;ele5};
end