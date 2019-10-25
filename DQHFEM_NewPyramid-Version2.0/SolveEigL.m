function [d,V,N]=SolveEigL(K,M,N,SE)

    [TN,~]=size(K);
    switch SE
        case 1
            N=TN;
            [V,Di]=eig(K,M);
            d=zeros(N,1);
            for j=1:N
                d(j,1)=Di(j,j);
            end

        case 2
            [V,Di]=eigs(K, M, N, 'sm'); 
            d=zeros(N,1);
            for j=1:N
                d(j,1)=Di(j,j);
            end
    end

end




