
function nq=OrderIntegral(number)
% 积分阶次：为了方便，不区分单元，取阶次最高的进行积分：（后期可优化）
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
    