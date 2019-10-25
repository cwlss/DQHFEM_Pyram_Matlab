function ElementData=CirPlate1(R,H,m,n)
% 金字塔+四面体=圆板
% 模型网格录入程序
ElementData=cell(8,1);
%%  Pyramid element 1-4
GX=cell(4,1);GY=GX;GZ=GX;
for ii=1:4
    X=cell(m-1,1);Y=X;Z=X;
    for i=1:m-1
        r=(m-i)*R/(m-1);
        theta=90*(ii-1):90/(m-i):90*ii;
        theta=theta';
        theta=theta/180*pi;
        xx=cos(theta)*r;
        yy=sin(theta)*r;
        zz=H*(m-i)/(m-1)*ones(size(xx));
        
        X{i}=xx;Y{i}=yy;Z{i}=zz;
    end
    node=cell(2,1);
    node{1}=[X{1}(1),Y{1}(1),H];
    node{2}=[X{1}(m),Y{1}(m),H];
    node{3}=[X{1}(m),Y{1}(m),0];
    node{4}=[X{1}(1),Y{1}(1),0];
    node{5}=[0,0,0];
    t=5;
    %边1（曲）
    for i=2:m-1
        t=t+1;
        node{t}=[X{1}(i),Y{1}(i),Z{1}(i)];
    end
    %     边2（直）
    for i=2:m-1
        t=t+1;
        node{t}=[X{1}(m),Y{1}(m),H*(m-i)/(m-1)];
    end
    %     边3（曲）
    for i=2:m-1
        t=t+1;
        node{t}=[X{1}(i),Y{1}(i),0];
    end
    %     边4（直）
    for i=2:m-1
        t=t+1;
        node{t}=[X{1}(1),Y{1}(1),H*(m-i)/(m-1)];
    end
    %     边5（曲）
    for i=2:m-1
        t=t+1;
        node{t}=[X{i}(1),Y{i}(1),Z{i}(1)];
    end
    %     边6（曲）
    for i=2:m-1
        t=t+1;
        k=length(X{i});
        node{t}=[X{i}(k),Y{i}(k),Z{i}(k)];
    end
    %     四边形面
    for i=2:m-1
        for j=2:m-1
            t=t+1;
            node{t}=[X{1}(j),Y{1}(j),H*(m-i)/(m-1)];
        end
    end
    %     三角形面2（曲）
    for i=2:m-2
        k=length(X{i});
        for j=2:k-1
            t=t+1;
            node{t}=[X{i}(j),Y{i}(j),Z{i}(j)];
        end
    end
    N=length(node);
    x1=zeros(N,1);y1=x1;z1=x1;
    for i=1:N
        x1(i)=node{i}(1);
        y1(i)=node{i}(2);
        z1(i)=node{i}(3);
    end
    GX{ii}=x1;
    GY{ii}=y1;
    GZ{ii}=z1;
end
for i=1:4
    ele=struct;
    ele.Type=4;
    number=struct;
    number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
    number.face={[n-1,n-1],n-2,n-2,n-2,n-2};
    number.H=n-2;
    ele.Number=number;
    number=struct;
    number.edge=[m-2,m-2,m-2,m-2,m-2,m-2,0,0];
    number.face={[m-2,m-2],m-3,0,0,0};
    number.H=0;
    ele.MapNumber=number;
    ele.MapRight={GX{i},GY{i},GZ{i}};
    ElementData{i}=ele;
end
%% Tetra element 5-8
GX=cell(4,1);GY=GX;GZ=GX;
for ii=1:4
    X=cell(m-1,1);Y=X;Z=X;
    for i=1:m-1
        r=(m-i)*R/(m-1);
        theta=90*(ii-1):90/(m-i):90*ii;
        theta=theta';
        theta=theta/180*pi;
        xx=cos(theta)*r;
        yy=sin(theta)*r;
        zz=H*(m-i)/(m-1)*ones(size(xx));
        X{i}=xx;Y{i}=yy;Z{i}=zz;
    end
    node=cell(2,1);
    node{1}=[0,0,H];
    node{2}=[X{1}(m),Y{1}(m),H];
    node{3}=[X{1}(1),Y{1}(1),H];
    node{4}=[0,0,0];
    t=4;
    for i=2:m-1
        t=t+1;
        node{t}=[X{1}(i),Y{1}(i),Z{1}(i)];
    end
    for i=m-1:-1:2
        t=t+1;
        k=length(X{i});
        node{t}=[X{i}(k),Y{i}(k),Z{i}(k)];
    end
    for i=m-1:-1:2
        t=t+1;
        node{t}=[X{i}(1),Y{i}(1),Z{i}(1)];
    end
    for i=2:m-2
        k=length(X{i});
        for j=2:k-1
            t=t+1;
            node{t}=[X{m+2-i-j}(j),Y{m+2-i-j}(j),Z{m+2-i-j}(j)];
        end
    end
    N=length(node);
    x1=zeros(N,1);y1=x1;z1=x1;
    for i=1:N
        x1(i)=node{i}(1);
        y1(i)=node{i}(2);
        z1(i)=node{i}(3);
    end
    GX{ii}=x1;
    GY{ii}=y1;
    GZ{ii}=z1;
end
for i=1:4
    ele=struct;
    ele.Type=3;
    number=struct;
    number.edge=[n-1,n-1,n-1,n-1,n-1,n-1];
    number.face={n-2,n-2,n-2,n-2};
    number.H=n-2;
    ele.Number=number;
    MapNumber.edge=[0,m-2,0,0,m-2,m-2];
    MapNumber.face={0,0,m-3,0};
    MapNumber.H=0;
    ele.MapNumber=MapNumber;
    ele.MapRight={GX{i},GY{i},GZ{i}};
    ElementData{i+4}=ele;
%     ElementData={ElementData{1}};
end
end
