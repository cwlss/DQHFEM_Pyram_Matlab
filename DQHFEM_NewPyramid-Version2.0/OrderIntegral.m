
function nq=OrderIntegral(number)
% ���ֽ״Σ�Ϊ�˷��㣬�����ֵ�Ԫ��ȡ�״���ߵĽ��л��֣������ڿ��Ż���
MAX=max(number.edge);
for i=1:length(number.face)
    Max=max(number.face{i});
    if MAX<Max
        MAX=Max;
    end
end
if MAX<max(number.H)
    MAX=max(number.H);
end
dq=1;
nq=[MAX+dq,MAX+dq,MAX+dq];
end
    