function bc=ApplyFaceBC(GC,GIH,FaceBC,Oriention)
% FaceBC={[a1,b1,c1,d1],...,[an,bn,cn,dn]}    an,bn,cn,dn为施加边界条件的平面的隐式方程Ax+By+Cz+D=0的系数A,B,C,D
%Oriention：
% 【1,1,1】——固支
% 【1,0,0】——X方向约束
% 【0,1,0】——Y方向约束
% 【0,0,1】——Z方向约束
%所以单向约束条件下，只给出了全局坐标的三个方向的约束。由于未建立局部坐标系，对于其余任意方向约束并不具备，待后期有时间建立局部坐标和全局坐标的转换关系。


N=length(GC{1});
bc=ones(3*N,1);
if nargin>2
    res=0.00000001;%点是否在给定面上的判定残差
    or=Oriention;
    for i=1:N
        if any(GIH==i)
        else
            for fi=1:length(FaceBC)
                A=FaceBC{fi}(1);
                B=FaceBC{fi}(2);
                C=FaceBC{fi}(3);
                D=FaceBC{fi}(4);
                if abs(A*GC{1}(i)+B*GC{2}(i)+C*GC{3}(i)+D)<res
                    if or(1)==1
                        bc(i)=0;
                    end
                    if or(2)==1
                        bc(i+N)=0;
                    end
                    if or(3)==1
                        bc(i+2*N)=0;
                    end
                end
            end
        end
    end
else
end
end
