function bc=ApplyFaceBC(GC,GIH,FaceBC,Oriention)
% FaceBC={[a1,b1,c1,d1],...,[an,bn,cn,dn]}    an,bn,cn,dnΪʩ�ӱ߽�������ƽ�����ʽ����Ax+By+Cz+D=0��ϵ��A,B,C,D
%Oriention��
% ��1,1,1��������֧
% ��1,0,0������X����Լ��
% ��0,1,0������Y����Լ��
% ��0,0,1������Z����Լ��
%���Ե���Լ�������£�ֻ������ȫ����������������Լ��������δ�����ֲ�����ϵ�������������ⷽ��Լ�������߱�����������ʱ�佨���ֲ������ȫ�������ת����ϵ��


N=length(GC{1});
bc=ones(3*N,1);
if nargin>2
    res=0.00000001;%���Ƿ��ڸ������ϵ��ж��в�
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
