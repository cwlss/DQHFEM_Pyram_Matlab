function ElementData=FiveTetra(a,b,h,n)
%% 模型录入程序
%%
n1=[0,0,0];
n2=[a,0,0];
n3=[a,b,0];
n4=[0,b,0];
n5=[0,0,h];
n6=[a,0,h];
n7=[a,b,h];
n8=[0,b,h];
% element 1 Ttra
ele1=struct;
ele1.Type=3;
ele1.MapRight={[n1(1);n2(1);n3(1);n6(1)],[n1(2);n2(2);n3(2);n6(2)],[n1(3);n2(3);n3(3);n6(3)]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele1.Number=number;


% element 2 Ttra
ele2=struct;
ele2.Type=3;
ele2.MapRight={[n5(1);n6(1);n1(1);n8(1)],[n5(2);n6(2);n1(2);n8(2)],[n5(3);n6(3);n1(3);n8(3)]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele2.Number=number;
% element 5 Ttra
ele3=struct;
ele3.Type=3;
ele3.MapRight={[n3(1);n4(1);n1(1);n8(1)],[n3(2);n4(2);n1(2);n8(2)],[n3(3);n4(3);n1(3);n8(3)]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele3.Number=number;

% element 4 Ttra
ele4=struct;
ele4.Type=3;
ele4.MapRight={[n6(1);n7(1);n3(1);n8(1)],[n6(2);n7(2);n3(2);n8(2)],[n6(3);n7(3);n3(3);n8(3)]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele4.Number=number;

% element 5 Ttra
ele5=struct;
ele5.Type=3;
ele5.MapRight={[n1(1);n3(1);n8(1);n6(1)],[n1(2);n3(2);n8(2);n6(2)],[n1(3);n3(3);n8(3);n6(3)]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele5.Number=number;




ElementData={ele1;ele2;ele3;ele4;ele5};
end

