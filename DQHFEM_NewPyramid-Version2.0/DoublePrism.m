function ElementData=DoublePrism(a,b,h,n)
%% 模型录入程序


% element2 Prism
ele2=struct;
ele2.Type=2;
ele2.MapRight={[0;a;0;0;a;0],[0;0;b;0;0;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele2.Number=number;

% element3 Prism
ele3=struct;
ele3.Type=2;
ele3.MapRight={[a;a;0;a;a;0],[0;b;b;0;b;b],[0;0;0;h;h;h]};
number=struct;
number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
number.face={n-2,[n-1,n-1],[n-1,n-1],[n-1,n-1],n-2};
number.H=[n-2,n-1];
ele3.Number=number;



ElementData={ele2;ele3};
end

    
    
    
    
    
    
    
