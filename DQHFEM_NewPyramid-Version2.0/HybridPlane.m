function ElementData=HybridPlane(a,b,h,n)
%% 模型录入程序
%%
%%
% element2 Prism
ele2=struct;
ele2.Type=2;
ele2.MapRight={[0;a/2;0;0;a/2;0],[0;0;b;0;0;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele2.Number=number;
% ele2.Vertex=ele2.MapRight;
MapNumber.edge=[0,0,0,0,0,0,0,0,0];
MapNumber.face={0,[0,0],[0,0],[0,0],0};
MapNumber.H=[0,0];
ele2.MapNumber=MapNumber;
%%
% element 3 Prism
ele3=struct;
ele3.Type=2;
ele3.MapRight={[a/2;a/2;0;a/2;a/2;0],[0;b;b;0;b;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele3.Number=number;
% ele3.Vertex=ele3.MapRight;
MapNumber.edge=[0,0,0,0,0,0,0,0,0];
MapNumber.face={0,[0,0],[0,0],[0,0],0};
MapNumber.H=[0,0];
ele3.MapNumber=MapNumber;
%%
% element 4 Pyramid
ele4=struct;
ele4.Type=4;
ele4.MapRight={[a/2;a/2;a/2;a/2;a],[0;0;b;b;0],[h;0;0;h;0]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],n-2,n-2,n-2,n-2};
number.H=n-2;
ele4.Number=number;
ele4.Vertex=ele4.MapRight;

MapNumber.edge=[0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],0,0,0,0};
MapNumber.H=0;
ele4.MapNumber=MapNumber;
%%
% element 5 Ttra
ele5=struct;
ele5.Type=3;
ele5.MapRight={[a;a;a/2;a/2],[0;0;b;0],[0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele5.Number=number;
% ele5.Vertex=ele5.MapRight;
MapNumber.edge=[0,0,0,0,0,0];
MapNumber.face={0,0,0,0};
MapNumber.H=0;
ele5.MapNumber=MapNumber;

% element 6 Pyramid
ele6=struct;
ele6.Type=4;
ele6.MapRight={[a/2;a;a;a/2;a],[b;b;b;b;0],[0;0;h;h;0]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],n-2,n-2,n-2,n-2};
number.H=n-2;
ele6.Number=number;
MapNumber.edge=[0,0,0,0,0,0,0,0];
MapNumber.face={[0,0],0,0,0,0};
MapNumber.H=0;
ele6.MapNumber=MapNumber;

% element 7 Ttra
ele7=struct;
ele7.Type=3;
ele7.MapRight={[a;a/2;a;a],[0;b;b;0],[h;h;h;0]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,n-2,n-2,n-2};
number.H=n-2;
ele7.Number=number;
MapNumber.edge=[0,0,0,0,0,0];
MapNumber.face={0,0,0,0};
MapNumber.H=0;
ele7.MapNumber=MapNumber;







ElementData={ele2;ele3;ele4;ele5;ele6;ele7};

end

    
    
    
    
    
    
    
