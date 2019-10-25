function ElementData=DoubleHex(a,b,h,n)
% 模型录入程序
%% element1 Hex
ele1=struct;
ele1.Type=1;
ele1.MapRight={[0;a/2;a/2;0;0;a/2;a/2;0],[0;0;b;b;0;0;b;b],[0;0;0;0;h;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1]};
number.H=[n-1,n-1,n-1];
ele1.Number=number;
ele1.Vertex=ele1.MapRight;

ele2=struct;
ele2.Type=1;
ele2.MapRight={[a/2;a;a;a/2;a/2;a;a;a/2],[0;0;b;b;0;0;b;b],[0;0;0;0;h;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1],[n-1,n-1]};
number.H=[n-1,n-1,n-1];
ele2.Number=number;
ele2.Vertex=ele2.MapRight;


% ElementData={ele1;ele2;ele3};
% ElementData={ele1;ele2;ele3;ele4;ele5};
% ElementData={ele4;ele5};
ElementData={ele1;ele2};
end

    
    
    
    
    
    
    
