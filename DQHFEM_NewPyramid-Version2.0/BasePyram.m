function [G,Gr,Gs,Gt]=BasePyramNN(x,y,z,number)
% 求取金字塔基的函数值和三个偏导数值
% number为结构体：储存金字塔单元边、面、体的阶次
%  可缺省调用 [G]=BasePyram(x,y,z,number)
%   此处显示详细说明
num=length(x);
ne=sum(number.edge);
nfq=number.face{1}(1)*number.face{1}(2);
nft=0;
for i=2:5
    tt=number.face{i};
    nt=tt*(tt+1)/2;
    nft=nft+nt;
end
H=number.H;
nh=0;
for i=1:H
    nh=nh+i*(i+1)/2;
end
NN=5+ne+nfq+nft+nh;
G=zeros(num,NN);

if nargout==1
    S=BaseVertex(x,y,z);
    G(:,1:5)=S;
    ni=5;
    for i=1:8
        nedgei=1:number.edge(i);
        PE=BaseEdge(x,y,z,number.edge(i),i);
        G(:,ni+nedgei)=PE;
        ni=ni+number.edge(i);
    end
    PF0=BaseFace(x,y,z,number.face{1},0);
    G(:,ni+1:ni+nfq)=PF0;
    ni=ni+nfq;
    
    for i=1:4
        fi=number.face{i+1};
        nfi=fi*(fi+1)/2;
        nfacei=1:nfi;
        PF=BaseFace(x,y,z,number.face{i+1},i);
        G(:,ni+nfacei)=PF;
        ni=ni+nfi;
    end
    PH=BaseH(x,y,z,number.H);
    G(:,ni+1:ni+nh)=PH;
else
    Gr=G;Gs=G;Gt=G;
    [S,dS]=BaseVertex(x,y,z);
    G(:,1:5)=S;
    Gr(:,1:5)=dS{1};
    Gs(:,1:5)=dS{2};
    Gt(:,1:5)=dS{3};
    ni=5;
    for i=1:8
        nedgei=1:number.edge(i);
        [PE,dPE]=BaseEdge(x,y,z,number.edge(i),i);
        G(:,ni+nedgei)=PE;
        Gr(:,ni+nedgei)=dPE{1};
        Gs(:,ni+nedgei)=dPE{2};
        Gt(:,ni+nedgei)=dPE{3};
        ni=ni+number.edge(i);
    end
    [PF0,dPF0]=BaseFace(x,y,z,number.face{1},0);
    G(:,ni+1:ni+nfq)=PF0;
    Gr(:,ni+1:ni+nfq)=dPF0{1};
    Gs(:,ni+1:ni+nfq)=dPF0{2};
    Gt(:,ni+1:ni+nfq)=dPF0{3};
    ni=ni+nfq;
    
    for i=1:4
        fi=number.face{i+1};
        nfi=fi*(fi+1)/2;
        nfacei=1:nfi;
        [PF,dPF]=BaseFace(x,y,z,number.face{i+1},i);
        G(:,ni+nfacei)=PF;
        Gr(:,ni+nfacei)=dPF{1};
        Gs(:,ni+nfacei)=dPF{2};
        Gt(:,ni+nfacei)=dPF{3};
        ni=ni+nfi;
    end
    [PH,dPH]=BaseH(x,y,z,number.H);
    G(:,ni+1:ni+nh)=PH;
    Gr(:,ni+1:ni+nh)=dPH{1};
    Gs(:,ni+1:ni+nh)=dPH{2};
    Gt(:,ni+1:ni+nh)=dPH{3};
end
end

function [S,dS] = BaseVertex(r,s,t)
nr=length(r);
S1=zeros(size(r));S2=S1;S3=S1;S4=S1;S5=S1;
dS1dr=S1;dS2dr=S1;dS3dr=S1;dS4dr=S1;dS5dr=S1;
dS1ds=S1;dS2ds=S1;dS3ds=S1;dS4ds=S1;dS5ds=S1;
dS1dt=S1;dS2dt=S1;dS3dt=S1;dS4dt=S1;dS5dt=S1;
for i=1:nr
    if t(i)==1
        S1(i)=0;
        S2(i)=0;
        S3(i)=0;
        S4(i)=0;
        S5(i)=1;
    else
        S1(i)=0.25*(1-r(i)-s(i)-t(i)+r(i)*s(i)/(1-t(i)));
        S2(i)=0.25*(1+r(i)-s(i)-t(i)-r(i)*s(i)/(1-t(i)));
        S3(i)=0.25*(1+r(i)+s(i)-t(i)+r(i)*s(i)/(1-t(i)));
        S4(i)=0.25*(1-r(i)+s(i)-t(i)-r(i)*s(i)/(1-t(i)));
        S5(i)=t(i);
    end
    if nargout>1
        if t(i)==1
            dS1dr(i)=0.25*(-1);   dS1ds(i)=0.25*(-1);   dS1dt(i)=0.25.*(-1);
            dS2dr(i)=0.25.*(1);   dS2ds(i)=0.25.*(-1);  dS2dt(i)=0.25.*(-1);
            dS3dr(i)=0.25.*(1);   dS3ds(i)=0.25.*(1);   dS3dt(i)=0.25.*(-1);
            dS4dr(i)=0.25.*(-1);  dS4ds(i)=0.25.*(1);   dS4dt(i)=0.25.*(-1);
            dS5dr(i)=0;           dS5ds(i)=0;           dS5dt(i)=1;
        else
            dS1dr(i)=0.25*(-1+s(i)/(1-t(i)));   dS1ds(i)=0.25*(-1+r(i)/(1-t(i)));   dS1dt(i)=0.25*(-1+r(i)*s(i)/(1-t(i)).^2);
            dS2dr(i)=0.25*(1-s(i)/(1-t(i)));    dS2ds(i)=0.25*(-1-r(i)/(1-t(i)));   dS2dt(i)=0.25*(-1-r(i)*s(i)/(1-t(i)).^2);
            dS3dr(i)=0.25*(1+s(i)/(1-t(i)));    dS3ds(i)=0.25*(1+r(i)/(1-t(i)));    dS3dt(i)=0.25*(-1+r(i)*s(i)/(1-t(i)).^2);
            dS4dr(i)=0.25*(-1-s(i)/(1-t(i)));   dS4ds(i)=0.25*(1-r(i)/(1-t(i)));    dS4dt(i)=0.25*(-1-r(i)*s(i)/(1-t(i)).^2);
            dS5dr(i)=0;                    dS5ds(i)=0;                    dS5dt(i)=1;
        end
    end
end
S=[S1,S2,S3,S4,S5];
if nargout>1
    Sr=[dS1dr,dS2dr,dS3dr,dS4dr,dS5dr];
    Ss=[dS1ds,dS2ds,dS3ds,dS4ds,dS5ds];
    St=[dS1dt,dS2dt,dS3dt,dS4dt,dS5dt];
    dS={Sr,Ss,St};
end

%%
end



function [P,dP] = BaseEdge(x,y,z,n,N)
%BaseEdge 此处显示有关此函数的摘要
%   此处显示详细说明
num=length(x);
PJ=zeros(num,n);PJx=PJ;PJy=PJ;PJz=PJ;
P=PJ;Px=PJ;Py=PJ;Pz=PJ;
if nargout==1
    Nv=BaseVertex(x,y,z);
    switch N
        case 1
            Poly=Nv(:,1).*Nv(:,2);
            PJ=BaseJacobiPyramR(n-1, 2, 2, x, z);
        case 2
            Poly=Nv(:,2).*Nv(:,3);
            PJ=BaseJacobiPyramS(n-1, 2, 2, y, z);
        case 3
            Poly=Nv(:,3).*Nv(:,4);
            PJ=BaseJacobiPyramR(n-1, 2, 2, x, z);
        case 4
            Poly=Nv(:,1).*Nv(:,4);
            PJ=BaseJacobiPyramS(n-1, 2, 2, y, z);
        case 5
            Poly=Nv(:,1).*Nv(:,5);
            for i=1:n
                PJ(:,i)=BaseJacobi(i-1,4,2, 2*z-1);
            end
        case 6
            Poly=Nv(:,2).*Nv(:,5);
            for i=1:n
                PJ(:,i)=BaseJacobi(i-1,4,2, 2*z-1);
            end
        case 7
            Poly=Nv(:,3).*Nv(:,5);
            for i=1:n
                PJ(:,i)=BaseJacobi(i-1,4,2, 2*z-1);
            end
        case 8
            Poly=Nv(:,4).*Nv(:,5);
            for i=1:n
                PJ(:,i)=BaseJacobi(i-1,4,2, 2*z-1);
            end
    end
    for i=1:n
        P(:,i)=Poly.*PJ(:,i);
    end
    
else  %the value of base and diff_base
    [Nv,dNv]=BaseVertex(x,y,z);
    dNvdr=dNv{1};dNvds=dNv{2};dNvdt=dNv{3};
    switch N
        case 1
            Poly=Nv(:,1).*Nv(:,2);
            Pdr=dNvdr(:,1).*Nv(:,2)+Nv(:,1).*dNvdr(:,2);
            Pds=dNvds(:,1).*Nv(:,2)+Nv(:,1).*dNvds(:,2);
            Pdt=dNvdt(:,1).*Nv(:,2)+Nv(:,1).*dNvdt(:,2);
            dPoly={Pdr,Pds,Pdt};
            [PJ,dPJ]=BaseJacobiPyramR(n-1, 2, 2, x, z);
            PJx=dPJ{1};PJy=dPJ{2};PJz=dPJ{3};
        case 2
            Poly=Nv(:,2).*Nv(:,3);
            Pdr=dNvdr(:,2).*Nv(:,3)+Nv(:,2).*dNvdr(:,3);
            Pds=dNvds(:,2).*Nv(:,3)+Nv(:,2).*dNvds(:,3);
            Pdt=dNvdt(:,2).*Nv(:,3)+Nv(:,2).*dNvdt(:,3);
            dPoly={Pdr,Pds,Pdt};
            [PJ,dPJ]=BaseJacobiPyramS(n-1, 2, 2, y, z);
            PJx=dPJ{1};PJy=dPJ{2};PJz=dPJ{3};
        case 3
            Poly=Nv(:,3).*Nv(:,4);
            Pdr=dNvdr(:,3).*Nv(:,4)+Nv(:,3).*dNvdr(:,4);
            Pds=dNvds(:,3).*Nv(:,4)+Nv(:,3).*dNvds(:,4);
            Pdt=dNvdt(:,3).*Nv(:,4)+Nv(:,3).*dNvdt(:,4);
            dPoly={Pdr,Pds,Pdt};
            [PJ,dPJ]=BaseJacobiPyramR(n-1, 2, 2, x, z);
            PJx=dPJ{1};PJy=dPJ{2};PJz=dPJ{3};
        case 4
            Poly=Nv(:,1).*Nv(:,4);
            Pdr=dNvdr(:,1).*Nv(:,4)+Nv(:,1).*dNvdr(:,4);
            Pds=dNvds(:,1).*Nv(:,4)+Nv(:,1).*dNvds(:,4);
            Pdt=dNvdt(:,1).*Nv(:,4)+Nv(:,1).*dNvdt(:,4);
            dPoly={Pdr,Pds,Pdt};
            [PJ,dPJ]=BaseJacobiPyramS(n-1, 2, 2, y, z);
            PJx=dPJ{1};PJy=dPJ{2};PJz=dPJ{3};
        case 5
            Poly=Nv(:,1).*Nv(:,5);
            Pdr=dNvdr(:,1).*Nv(:,5)+Nv(:,1).*dNvdr(:,5);
            Pds=dNvds(:,1).*Nv(:,5)+Nv(:,1).*dNvds(:,5);
            Pdt=dNvdt(:,1).*Nv(:,5)+Nv(:,1).*dNvdt(:,5);
            dPoly={Pdr,Pds,Pdt};
            for i=1:n
                [PJ(:,i),dPJ]=BaseJacobi(i-1,4,2, 2*z-1);
                PJx(:,i)=0*dPJ;PJy(:,i)=0*dPJ;PJz(:,i)=2*dPJ;
            end
        case 6
            Poly=Nv(:,5).*Nv(:,2);
            Pdr=dNvdr(:,5).*Nv(:,2)+Nv(:,5).*dNvdr(:,2);
            Pds=dNvds(:,5).*Nv(:,2)+Nv(:,5).*dNvds(:,2);
            Pdt=dNvdt(:,5).*Nv(:,2)+Nv(:,5).*dNvdt(:,2);
            dPoly={Pdr,Pds,Pdt};
            for i=1:n
                [PJ(:,i),dPJ]=BaseJacobi(i-1,4,2, 2*z-1);
                PJx(:,i)=0*dPJ;PJy(:,i)=0*dPJ;PJz(:,i)=2*dPJ;
            end
        case 7
            Poly=Nv(:,3).*Nv(:,5);
            Pdr=dNvdr(:,3).*Nv(:,5)+Nv(:,3).*dNvdr(:,5);
            Pds=dNvds(:,3).*Nv(:,5)+Nv(:,3).*dNvds(:,5);
            Pdt=dNvdt(:,3).*Nv(:,5)+Nv(:,3).*dNvdt(:,5);
            dPoly={Pdr,Pds,Pdt};
            for i=1:n
                [PJ(:,i),dPJ]=BaseJacobi(i-1,4,2, 2*z-1);
                PJx(:,i)=0*dPJ;PJy(:,i)=0*dPJ;PJz(:,i)=2*dPJ;
            end
        case 8
            Poly=Nv(:,4).*Nv(:,5);
            Pdr=dNvdr(:,4).*Nv(:,5)+Nv(:,4).*dNvdr(:,5);
            Pds=dNvds(:,4).*Nv(:,5)+Nv(:,4).*dNvds(:,5);
            Pdt=dNvdt(:,4).*Nv(:,5)+Nv(:,4).*dNvdt(:,5);
            dPoly={Pdr,Pds,Pdt};
            for i=1:n
                [PJ(:,i),dPJ]=BaseJacobi(i-1,4,2, 2*z-1);
                PJx(:,i)=0*dPJ;PJy(:,i)=0*dPJ;PJz(:,i)=2*dPJ;
            end
    end
    for i=1:n
        P(:,i)=Poly.*PJ(:,i);
        Px(:,i)=dPoly{1}.*PJ(:,i)+Poly.*PJx(:,i);
        Py(:,i)=dPoly{2}.*PJ(:,i)+Poly.*PJy(:,i);
        Pz(:,i)=dPoly{3}.*PJ(:,i)+Poly.*PJz(:,i);
    end
    dP={Px,Py,Pz};
end
end
function [P,dP] = BaseFace(x,y,z,n,N)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
if nargout==1
    switch N
        case 0
            if length(n)~=2
                error('n的输入长度有误');
            end
            P=BaseFace0(x,y,z,n);
        case 1
 
            P=BaseFace13( x,y,z,n,1 );
        case 2
            P=BaseFace24( x,y,z,n,2 );
        case 3
            P=BaseFace13( x,y,z,n,3 );
        case 4
            P=BaseFace24( x,y,z,n,4 );
    end
else
    switch N
        case 0
            if length(n)~=2
                error('n的输入长度有误');
            end
            [P,dP]=BaseFace0(x,y,z,n);
        case 1
            [P,dP]=BaseFace13( x,y,z,n,1 );
        case 2
            [P,dP]=BaseFace24( x,y,z,n,2);
        case 3
            [P,dP]=BaseFace13( x,y,z,n,3 );
        case 4
            [P,dP]=BaseFace24( x,y,z,n,4 );
    end
end
end
function [P,dP]=BaseFace13(x,y,z,n,Type)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
num=length(x);
nft=n*(n+1)/2;
P=zeros(num,nft);
if nargout==1
    %
    Nv=BaseVertex(x,y,z);
    if Type==1
        Poly=Nv(:,1).*Nv(:,2).*Nv(:,5);
    elseif Type==3
        Poly=Nv(:,3).*Nv(:,4).*Nv(:,5);
    else
        error('error');
    end
    %
    ij=1;
    Pi=BaseJacobiPyramR(n-1, 2, 2, x, z);
    for i=1:n
        for j=1:n-i+1
            Pj=BaseJacobi(j-1, 2*(i-1)+6, 2, 2*z-1);
            P(:,ij)=Poly.*Pi(:,i).*Pj;
            ij=ij+1;
        end
    end
else
    Px=zeros(num,nft);Py=Px;Pz=Px;
    %
    [Nv,dNv]=BaseVertex(x,y,z);
    dNvdr=dNv{1};dNvds=dNv{2};dNvdt=dNv{3};
    switch Type
        case 1
        Poly=Nv(:,1).*Nv(:,2).*Nv(:,5);
        Pdr=dNvdr(:,1).*Nv(:,2).*Nv(:,5)+Nv(:,1).*dNvdr(:,2).*Nv(:,5)+Nv(:,1).*Nv(:,2).*dNvdr(:,5);
        Pds=dNvds(:,1).*Nv(:,2).*Nv(:,5)+Nv(:,1).*dNvds(:,2).*Nv(:,5)+Nv(:,1).*Nv(:,2).*dNvds(:,5);
        Pdt=dNvdt(:,1).*Nv(:,2).*Nv(:,5)+Nv(:,1).*dNvdt(:,2).*Nv(:,5)+Nv(:,1).*Nv(:,2).*dNvdt(:,5);
        dPoly={Pdr,Pds,Pdt};
        case 3
        Poly=Nv(:,3).*Nv(:,4).*Nv(:,5);
        Pdr=dNvdr(:,3).*Nv(:,4).*Nv(:,5)+Nv(:,3).*dNvdr(:,4).*Nv(:,5)+Nv(:,3).*Nv(:,4).*dNvdr(:,5);
        Pds=dNvds(:,3).*Nv(:,4).*Nv(:,5)+Nv(:,3).*dNvds(:,4).*Nv(:,5)+Nv(:,3).*Nv(:,4).*dNvds(:,5);
        Pdt=dNvdt(:,3).*Nv(:,4).*Nv(:,5)+Nv(:,3).*dNvdt(:,4).*Nv(:,5)+Nv(:,3).*Nv(:,4).*dNvdt(:,5);
        dPoly={Pdr,Pds,Pdt};
    end
% % 
    ij=1;
    [Pi,dPi]=BaseJacobiPyramR(n-1, 2, 2, x, z);
    for i=1:n
        for j=1:n-i+1
            [Pj,dPj]=BaseJacobi(j-1, 2*(i-1)+6, 2, 2*z-1);
            P(:,ij)=Poly.*Pi(:,i).*Pj;
            Px(:,ij)=dPoly{1}.*Pi(:,i).*Pj+Poly.*dPi{1}(:,i).*Pj+Poly.*Pi(:,i).*0.*dPj;
            Py(:,ij)=dPoly{2}.*Pi(:,i).*Pj+Poly.*dPi{2}(:,i).*Pj+Poly.*Pi(:,i).*0.*dPj;
            Pz(:,ij)=dPoly{3}.*Pi(:,i).*Pj+Poly.*dPi{3}(:,i).*Pj+Poly.*Pi(:,i).*2.*dPj;
            ij=ij+1;
        end
    end
    dP={Px,Py,Pz};
end
end
function [P,dP]=BaseFace24(x,y,z,n,Type)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
num=length(x);
nft=n*(n+1)/2;
P=zeros(num,nft);
if nargout==1
    %
    Nv=BaseVertex(x,y,z);
    if Type==2
        Poly=Nv(:,2).*Nv(:,3).*Nv(:,5);
    elseif Type==4
        Poly=Nv(:,1).*Nv(:,4).*Nv(:,5);
    else
        error('error');
    end
    %
    ij=1;
    Pi=BaseJacobiPyramS(n-1, 2, 2, y, z);
    for i=1:n
        for j=1:n-i+1
            Pj=BaseJacobi(j-1, 2*(i-1)+6, 2, 2*z-1);
            P(:,ij)=Poly.*Pi(:,i).*Pj;
            ij=ij+1;
        end
    end
else
    Px=zeros(num,nft);Py=Px;Pz=Px;
    %
    [Nv,dNv]=BaseVertex(x,y,z);
    dNvdr=dNv{1};dNvds=dNv{2};dNvdt=dNv{3};
    switch Type
        case 2
        Poly=Nv(:,3).*Nv(:,2).*Nv(:,5);
        Pdr=dNvdr(:,3).*Nv(:,2).*Nv(:,5)+Nv(:,3).*dNvdr(:,2).*Nv(:,5)+Nv(:,3).*Nv(:,2).*dNvdr(:,5);
        Pds=dNvds(:,3).*Nv(:,2).*Nv(:,5)+Nv(:,3).*dNvds(:,2).*Nv(:,5)+Nv(:,3).*Nv(:,2).*dNvds(:,5);
        Pdt=dNvdt(:,3).*Nv(:,2).*Nv(:,5)+Nv(:,3).*dNvdt(:,2).*Nv(:,5)+Nv(:,3).*Nv(:,2).*dNvdt(:,5);
        dPoly={Pdr,Pds,Pdt};
        case 4
        Poly=Nv(:,1).*Nv(:,4).*Nv(:,5);
        Pdr=dNvdr(:,1).*Nv(:,4).*Nv(:,5)+Nv(:,1).*dNvdr(:,4).*Nv(:,5)+Nv(:,1).*Nv(:,4).*dNvdr(:,5);
        Pds=dNvds(:,1).*Nv(:,4).*Nv(:,5)+Nv(:,1).*dNvds(:,4).*Nv(:,5)+Nv(:,1).*Nv(:,4).*dNvds(:,5);
        Pdt=dNvdt(:,1).*Nv(:,4).*Nv(:,5)+Nv(:,1).*dNvdt(:,4).*Nv(:,5)+Nv(:,1).*Nv(:,4).*dNvdt(:,5);
        dPoly={Pdr,Pds,Pdt};
    end
% % 
    ij=1;
    [Pi,dPi]=BaseJacobiPyramS(n-1, 2, 2, y, z);
    for i=1:n
        for j=1:n-i+1
            [Pj,dPj]=BaseJacobi(j-1, 2*(i-1)+6, 2, 2*z-1);
            P(:,ij)=Poly.*Pi(:,i).*Pj;
            Px(:,ij)=dPoly{1}.*Pi(:,i).*Pj+Poly.*dPi{1}(:,i).*Pj+Poly.*Pi(:,i).*0.*dPj;
            Py(:,ij)=dPoly{2}.*Pi(:,i).*Pj+Poly.*dPi{2}(:,i).*Pj+Poly.*Pi(:,i).*0.*dPj;
            Pz(:,ij)=dPoly{3}.*Pi(:,i).*Pj+Poly.*dPi{3}(:,i).*Pj+Poly.*Pi(:,i).*2.*dPj;
            ij=ij+1;
        end
    end
    dP={Px,Py,Pz};
end
end
function [P,dP]=BaseFace0(x,y,z,n)
%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
nfq=n(1)*n(2);
num=length(x);
if nargout==1
    P=zeros(num,nfq);
%     
    Nv=BaseVertex(x,y,z);
    Poly=Nv(:,1).*Nv(:,3);
%     
    ij=1;
    Pj=BaseJacobiPyramS(n(2)-1, 2, 2, y, z);
    Pi=BaseJacobiPyramR(n(1)-1, 2, 2, x, z);
    for j=1:n(2)
        for i=1:n(1)
            P(:,ij)=Poly.*Pi(:,i).*Pj(:,j);
            ij=ij+1;
        end
    end
else
    P=zeros(num,nfq);Px=P;Py=P;Pz=P;
    
    [Nv,dNv]=BaseVertex(x,y,z);
    dNvdr=dNv{1};dNvds=dNv{2};dNvdt=dNv{3};
    Poly=Nv(:,1).*Nv(:,3);
    Pdr=dNvdr(:,1).*Nv(:,3)+Nv(:,1).*dNvdr(:,3);
    Pds=dNvds(:,1).*Nv(:,3)+Nv(:,1).*dNvds(:,3);
    Pdt=dNvdt(:,1).*Nv(:,3)+Nv(:,1).*dNvdt(:,3);
    dPoly={Pdr,Pds,Pdt};
    
    ij=1;
    [Pj,dPj]=BaseJacobiPyramS(n(2)-1, 2, 2, y, z);
    [Pi,dPi]=BaseJacobiPyramR(n(1)-1, 2, 2, x, z);
    for j=1:n(2)
        for i=1:n(1)
            P(:,ij)=Poly.*Pi(:,i).*Pj(:,j);
            Px(:,ij)=dPoly{1}.*Pi(:,i).*Pj(:,j)+Poly.*dPi{1}(:,i).*Pj(:,j)+Poly.*Pi(:,i).*dPj{1}(:,j);
            Py(:,ij)=dPoly{2}.*Pi(:,i).*Pj(:,j)+Poly.*dPi{2}(:,i).*Pj(:,j)+Poly.*Pi(:,i).*dPj{2}(:,j);
            Pz(:,ij)=dPoly{3}.*Pi(:,i).*Pj(:,j)+Poly.*dPi{3}(:,i).*Pj(:,j)+Poly.*Pi(:,i).*dPj{3}(:,j);
            ij=ij+1;
        end
    end
    dP={Px,Py,Pz};
end
end
function [P,dP]=BaseH(x,y,z,n)

%UNTITLED2 此处显示有关此函数的摘要
%   此处显示详细说明
% if length(n)~=3
%     error('n输入有误')
% end
% num=length(x);
% nh=n(1)*n(2)*n(3);
%
% if nargout==1
%     P=zeros(num,nh);
%     pe=[1,1,1,1,1];
%     Poly=BasePolyPart(x,y,z,pe);
%     ijk=1;
%     Pi=BaseJacobiPyramR(n(1)-1, 2, 2, x, z);
%     Pj=BaseJacobiPyramS(n(2)-1, 2, 2, y, z);
%     for i=1:n(1)
%         for j=1:n(2)
%             for k=1:n(3)
%                 Pk=BaseJacobi(k-1, 2*(i-1+j-1+5), 2, 2*z-1);
%                 P(:,ijk)=Poly.*Pi(:,i).*Pj(:,j).*Pk;
%                 ijk=ijk+1;
%             end
%         end
%     end
% else
%     P=zeros(num,nh);Px=P;Py=P;Pz=P;
%     pe=[1,1,1,1,1];
%     [Poly,dPoly]=BasePolyPart(x,y,z,pe);
%     ijk=1;
%     [Pi,dPi]=BaseJacobiPyramR(n(1)-1, 2, 2, x, z);
%     [Pj,dPj]=BaseJacobiPyramS(n(2)-1, 2, 2, y, z);
%     for i=1:n(1)
%         for j=1:n(2)
%             for k=1:n(3)
%                 [Pk,dPk]=BaseJacobi(k-1, 2*(i-1+j-1+5), 2, 2*z-1);
%                 P(:,ijk)=Poly.*Pi(:,i).*Pj(:,j).*Pk;
%                 Px(:,ijk)=dPoly{1}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{1}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{1}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*0.*dPk;
%                 Py(:,ijk)=dPoly{2}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{2}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{2}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*0.*dPk;
%                 Pz(:,ijk)=dPoly{3}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{3}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{3}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*2.*dPk;
%                 ijk=ijk+1;
%             end
%         end
%     end
%     dP={Px,Py,Pz};
% end
num=length(x);
nh=0;
for i=1:n
    nh=nh+i*(i+1)/2;
end
if nargout==1
    Nv=zeros(num,5);
    for i=1:num
        if z(i)==1
            Nv(i,1)=0;
            Nv(i,2)=0;
            Nv(i,3)=0;
            Nv(i,4)=0;
            Nv(i,5)=1;
        else
            Nv(i,1)=0.25.*(1-x(i)-y(i)-z(i)+x(i).*y(i)./(1-z(i)));
            Nv(i,2)=0.25.*(1+x(i)-y(i)-z(i)-x(i).*y(i)./(1-z(i)));
            Nv(i,3)=0.25.*(1+x(i)+y(i)-z(i)+x(i).*y(i)./(1-z(i)));
            Nv(i,4)=0.25.*(1-x(i)+y(i)-z(i)-x(i).*y(i)./(1-z(i)));
            Nv(i,5)=z(i);
        end
    end
    P=zeros(num,nh);
    %     pe=[1,1,1,1,1];
    %     Poly=BasePolyPart(x,y,z,pe);
    Poly=Nv(:,1).*Nv(:,3).*Nv(:,5);
    ijk=1;
    Pi=BaseJacobiPyramR(n-1, 2, 2, x, z);
    Pj=BaseJacobiPyramS(n-1, 2, 2, y, z);
    for i=1:n
        for j=1:n-i+1
            for k=1:n-i-j+2
                Pk=BaseJacobi(k-1, 2*(i-1+j-1+3), 2, 2*z-1);
                P(:,ijk)=Poly.*Pi(:,i).*Pj(:,j).*Pk;
                ijk=ijk+1;
            end
        end
    end
else
     Nv=zeros(num,5);dNvdr=zeros(num,5);dNvds=zeros(num,5);dNvdt=zeros(num,5);
    for i=1:num
        if z(i)==1
            dNvdr(i,1)=0.25*(-1);   dNvds(i,1)=0.25*(-1);  dNvdt(i,1)=0.25.*(-1);Nv(i,1)=0;
            dNvdr(i,2)=0.25.*(1);   dNvds(i,2)=0.25.*(-1); dNvdt(i,2)=0.25.*(-1);Nv(i,2)=0;
            dNvdr(i,3)=0.25.*(1);   dNvds(i,3)=0.25.*(1);  dNvdt(i,3)=0.25.*(-1);Nv(i,3)=0;
            dNvdr(i,4)=0.25.*(-1);  dNvds(i,4)=0.25.*(1);  dNvdt(i,4)=0.25.*(-1);Nv(i,4)=0;
            dNvdr(i,5)=0;           dNvds(i,5)=0;          dNvdt(i,5)=1;         Nv(i,5)=1;
        else
            dNvdr(i,1)=0.25*(-1+y(i)./(1-z(i)));
            dNvds(i,1)=0.25*(-1+x(i)./(1-z(i)));
            dNvdt(i,1)=0.25*(-1+x(i)*y(i)./(1-z(i)).^2);
            Nv(i,1)=0.25.*(1-x(i)-y(i)-z(i)+x(i).*y(i)./(1-z(i)));
            
            dNvdr(i,2)=0.25*(1-y(i)./(1-z(i)));
            dNvds(i,2)=0.25*(-1-x(i)./(1-z(i)));
            dNvdt(i,2)=0.25*(-1-x(i)*y(i)./(1-z(i)).^2);
            Nv(i,2)=0.25.*(1+x(i)-y(i)-z(i)-x(i).*y(i)./(1-z(i)));
            
            dNvdr(i,3)=0.25*(1+y(i)./(1-z(i)));
            dNvds(i,3)=0.25*(1+x(i)./(1-z(i)));
            dNvdt(i,3)=0.25*(-1+x(i)*y(i)./(1-z(i)).^2);
            Nv(i,3)=0.25.*(1+x(i)+y(i)-z(i)+x(i).*y(i)./(1-z(i)));
            
            dNvdr(i,4)=0.25*(-1-y(i)./(1-z(i)));
            dNvds(i,4)=0.25*(1-x(i)./(1-z(i)));
            dNvdt(i,4)=0.25*(-1-x(i)*y(i)./(1-z(i)).^2);
            Nv(i,4)=0.25.*(1-x(i)+y(i)-z(i)-x(i).*y(i)./(1-z(i)));
            
            dNvdr(i,5)=0;
            dNvds(i,5)=0;
            dNvdt(i,5)=1;
            Nv(i,5)=z(i);
        end
    end
    P=zeros(num,nh);Px=P;Py=P;Pz=P;
%     pe=[1,1,1,1,1];
%     [Poly,dPoly]=BasePolyPart(x,y,z,pe);
     Poly=Nv(:,1).*Nv(:,3).*Nv(:,5);
     dPoly{1}=dNvdr(:,1).*Nv(:,3).*Nv(:,5)+Nv(:,1).*dNvdr(:,3).*Nv(:,5)+Nv(:,1).*Nv(:,3).*dNvdr(:,5);
     dPoly{2}=dNvds(:,1).*Nv(:,3).*Nv(:,5)+Nv(:,1).*dNvds(:,3).*Nv(:,5)+Nv(:,1).*Nv(:,3).*dNvds(:,5);
     dPoly{3}=dNvdt(:,1).*Nv(:,3).*Nv(:,5)+Nv(:,1).*dNvdt(:,3).*Nv(:,5)+Nv(:,1).*Nv(:,3).*dNvdt(:,5);

   
     ijk=1;
    [Pi,dPi]=BaseJacobiPyramR(n-1, 2, 2, x, z);
    [Pj,dPj]=BaseJacobiPyramS(n-1, 2, 2, y, z);
    for i=1:n
        for j=1:n-i+1
            for k=1:n-i-j+2
                [Pk,dPk]=BaseJacobi(k-1, 2*(i-1+j-1+3), 2, 2*z-1);
                P(:,ijk)=Poly.*Pi(:,i).*Pj(:,j).*Pk;
                Px(:,ijk)=dPoly{1}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{1}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{1}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*0.*dPk;
                Py(:,ijk)=dPoly{2}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{2}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{2}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*0.*dPk;
                Pz(:,ijk)=dPoly{3}.*Pi(:,i).*Pj(:,j).*Pk+Poly.*dPi{3}(:,i).*Pj(:,j).*Pk+Poly.*Pi(:,i).*dPj{3}(:,j).*Pk+Poly.*Pi(:,i).*Pj(:,j).*2.*dPk;
                ijk=ijk+1;
            end
        end
    end
    dP={Px,Py,Pz};
end
end

function [Pi, jac]=BaseJacobiPyramR(N, a, b, x, z)

% Get the Jacobi polynomials and their derivatives on a pyramid
%
%  Calling sequence:
%
%    Pi=JacobiPyrami(N, a, b, x, z)
%
%    [Pi, jac]=JacobiPyrami(N, a, b, x, z)
%
%  Input :
%
%    N = The index of basis
%
%    a, b - Parameters of Jacobi polynomials
%
%    x, z - A vector of evaluation points on a unit pyramid
%
%  Output :
%
%    Pi - The Jacobi polynomials on a unit pyramid
%
%    jac - The first derivatives of Jacobi polynomials on a unit pyramid
%

% Coodinate transformation
s=x; t=1-z;

% Prepare matrices
n=length(s);
Pi=zeros(n, N+1);
if nargout>=2
    jac{1}=Pi; jac{2}=Pi; jac{3}=Pi;
end

% Evaluate
if N==0
    Pi(:,1)=ones(n, 1);
elseif N==1
    P=t*(a-b)/2+s*(a+b+2)/2;
    Pi(:,1)=ones(n, 1);
    Pi(:,2)=P;
    if nargout>=2
        dPs=(a+b+2)/2;
        dPt=(a-b)/2;
        jac{1}(:,2)=dPs*Pi(:,1);
        jac{2}(:,2)=0*Pi(:,1);
        jac{3}(:,2)=(-dPt)*Pi(:,1);
    end
else
    P1=1; P2=t*(a-b)/2+s*(a+b+2)/2;
    Pi(:,1)=ones(n, 1); Pi(:,2)=P2;
    if nargout>=2
        dPs1=0; dPs2=(a+b+2)/2;
        dPt1=0; dPt2=(a-b)/2;
        jac{1}(:,2)=dPs2*Pi(:,1);
        jac{2}(:,2)=0*Pi(:,1);
        jac{3}(:,2)=(-dPt2)*Pi(:,1);
    end
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b);
        a2=(2*n+a+b+1)*(a^2-b^2);
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
        a4=2*(n+a)*(n+b)*(2*n+a+b+2);
        P=((a2*t+a3*s).*P2-a4*t.^2.*P1)/a1;
        Pi(:,n+2)=P;
        if nargout>=2
            dPs=((a2*t+a3*s).*dPs2+a3*P2-a4*t.^2.*dPs1)/a1;
            dPt=((a2*t+a3*s).*dPt2+a2*P2-a4*t.^2.*dPt1-2*a4*t.*P1)/a1;
            jac{1}(:,n+2)=dPs;
            jac{2}(:,n+2)=0*dPs;
            jac{3}(:,n+2)=-dPt;
        end
        
        P1=P2; P2=P;
        if nargout>=2
            dPs1=dPs2; dPs2=dPs;
            dPt1=dPt2; dPt2=dPt;
        end
    end
end

%% demo
% % The node number for plot (N) and parameters of Jacobi polynomials (a, b)
% N=13; a=3; b=2; i=3;
%
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
%
% % Get nodes
% [x, y, z]=PyramLobatto(N);
% tri=pyramdelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)'];
%
% % Get the basis on the unti tetrahedron
% [Pi, jac]=JacobiPyrami(i, a, b, x, z);
% dPix=jac{1}; dPiy=jac{2}; dPiz=jac{3};
%
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% dPixi=Dx*Pi;
% dPiyi=Dy*Pi;
% dPizi=Dz*Pi;
%
% % Draw basis
% k=3;
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, Pi(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
%
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, dPix(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
%
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, dPixi(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [Pi, jac]=BaseJacobiPyramS(N, a, b, y, z)

% Get the Jacobi polynomials and their derivatives on a pyramid
%
%  Calling sequence:
%
%    Pi=JacobiPyramj(N, a, b, x, z)
%
%    [Pi, jac]=JacobiPyramj(N, a, b, x, z)
%
%  Input :
%
%    N = The index of basis
%
%    a, b - Parameters of Jacobi polynomials
%
%    x, z - A vector of evaluation points on a unit pyramid
%
%  Output :
%
%    Pi - The Jacobi polynomials on a unit pyramid
%
%    jac - The first derivatives of Jacobi polynomials on a unit pyramid
%

% Coodinate transformation
s=y; t=1-z;

% Prepare matrices
n=length(s);
Pi=zeros(n, N+1);
if nargout>=2
    jac{1}=Pi; jac{2}=Pi; jac{3}=Pi;
end

% Evaluate
if N==0
    Pi(:,1)=ones(n, 1);
elseif N==1
    P=t*(a-b)/2+s*(a+b+2)/2;
    Pi(:,1)=ones(n, 1);
    Pi(:,2)=P;
    if nargout>=2
        dPs=(a+b+2)/2;
        dPt=(a-b)/2;
        jac{1}(:,2)=0*Pi(:,1);
        jac{2}(:,2)=dPs*Pi(:,1);
        jac{3}(:,2)=(-dPt)*Pi(:,1);
    end
else
    P1=1; P2=t*(a-b)/2+s*(a+b+2)/2;
    Pi(:,1)=ones(n, 1); Pi(:,2)=P2;
    if nargout>=2
        dPs1=0; dPs2=(a+b+2)/2;
        dPt1=0; dPt2=(a-b)/2;
        jac{1}(:,2)=0*Pi(:,1);
        jac{2}(:,2)=dPs2*Pi(:,1);
        jac{3}(:,2)=(-dPt2)*Pi(:,1);
    end
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b);
        a2=(2*n+a+b+1)*(a^2-b^2);
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
        a4=2*(n+a)*(n+b)*(2*n+a+b+2);
        P=((a2*t+a3*s).*P2-a4*t.^2.*P1)/a1;
        Pi(:,n+2)=P;
        if nargout>=2
            dPs=((a2*t+a3*s).*dPs2+a3*P2-a4*t.^2.*dPs1)/a1;
            dPt=((a2*t+a3*s).*dPt2+a2*P2-a4*t.^2.*dPt1-2*a4*t.*P1)/a1;
            jac{1}(:,n+2)=0*dPs;
            jac{2}(:,n+2)=dPs;
            jac{3}(:,n+2)=-dPt;
        end
        
        P1=P2; P2=P;
        if nargout>=2
            dPs1=dPs2; dPs2=dPs;
            dPt1=dPt2; dPt2=dPt;
        end
    end
end


%% demo
% % The node number for plot (N) and parameters of Jacobi polynomials (a, b)
% N=13; a=3; b=2; i=3;
%
% % The number of points covered by RBF
% nup=3*N; % this parameter can be optimized ...
%
% % Get nodes
% [x, y, z]=PyramLobatto(N);
% tri=pyramdelaunay(N);
% pnts=[x(:)'; y(:)'; z(:)'];
%
% % Get the basis on the unti tetrahedron
% [Pi, jac]=JacobiPyramj(i, a, b, y, z);
% dPix=jac{1}; dPiy=jac{2}; dPiz=jac{3};
%
% % Radial basis 3D
% [ARBF, Dx, Dy, Dz]=rbfbasis3d(pnts, pnts, nup);
% dPixi=Dx*Pi;
% dPiyi=Dy*Pi;
% dPizi=Dz*Pi;
%
% % Draw basis
% k=3;
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, Pi(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
%
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, dPiy(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('Exact');
%
% figure; hold on;
% for i=1:N+3
%     trisurf(tri{i}, x, y, z, dPiyi(:,k));
% end
% shading interp;
% colorbar;
% view(3); alpha(0.5);
% axis equal;
% title('RBF');
end
function [ P,dP] =BasePolyPart(x,y,z,pe)
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
A=pe(1);B=pe(2);C=pe(3);D=pe(4);E=pe(5);
a=z.^A;
b=(1-x-z).^B;
c=(1-y-z).^C;
d=(1+x-z).^D;
e=(1+y-z).^E;
if nargout==1
    P=a.*b.*c.*d.*e;
else
    if A==0
        da=0*z;
    else
        da=A*z.^(A-1);
    end
    if B==0
        db=0*z;
    else
        db=B*(1-x-z).^(B-1);
    end
    if C==0
        dc=0*z;
    else
        dc=C*(1-y-z).^(C-1);
    end
    if D==0
        dd=0*z;
    else
        dd=D*(1+x-z).^(D-1);
    end
    if E==0
        de=0*z;
    else
        de=E*(1+y-z).^(E-1);
    end
    P=a.*b.*c.*d.*e;
    Px=a.*c.*e.*(-db.*d+dd.*b);
    Py=a.*b.*d.*(-dc.*e+de.*c);
    Pz=da.*b.*c.*d.*e-a.*db.*c.*d.*e-a.*b.*dc.*d.*e-a.*b.*c.*dd.*e-a.*b.*c.*d.*de;
    dP={Px,Py,Pz};
end
end
function [y, dy, ddy, dddy]=BaseJacobi(N, a, b, x)

% Get the Jacobi polynomials and their derivatives
%
%  Calling sequence:
%
%    y=JacobiRecDer(N, a, b, x)
%
%    [y, dy]=JacobiRecDer(N, a, b, x)
%
%    [y, dy, ddy]=JacobiRecDer(N, a, b, x)
%
%    [y, dy, ddy, dddy]=JacobiRecDer(N, a, b, x)
%
%  Input :
%
%    N = the number of basis
%    a, b - parameters of Jacobi polynomials
%    x - a vector of nodes
%
%  Output :
%
%    y - the Jacobi polynomials
%    dy, ddy, dddy - the first, second and third
%          order derivatives of Jacobi polynomials
%

if N<0
    y=zeros(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
elseif N==0
    y=ones(size(x)); dy=0*y; ddy=0*y; dddy=0*y;
elseif N==1
    y=0.5*(a-b+(a+b+2)*x); dy=0.5*(a+b+2); ddy=0*y; dddy=0*y;
else
    y1=1; y2=0.5*(a-b+(a+b+2)*x);
    dy1=0; dy2=0.5*(a+b+2);
    ddy1=0; ddy2=0;
    dddy1=0; dddy2=0;
    for n=1:N-1
        a1=2*(n+1)*(n+a+b+1)*(2*n+a+b);
        a2=(2*n+a+b+1)*(a^2-b^2);
        a3=(2*n+a+b)*(2*n+a+b+1)*(2*n+a+b+2);
        a4=2*(n+a)*(n+b)*(2*n+a+b+2);
        y=((a2+a3*x).*y2-a4*y1)/a1;
        if nargout>1
            dy=((a2+a3*x).*dy2+a3*y2-a4*dy1)/a1;
        end
        if nargout>2
            ddy=((a2+a3*x).*ddy2+2*a3*dy2-a4*ddy1)/a1;
        end
        if nargout>3
            dddy=((a2+a3*x).*dddy2+3*a3*ddy2-a4*dddy1)/a1;
        end
        
        y1=y2; y2=y;
        if nargout>1
            dy1=dy2; dy2=dy;
        end
        if nargout>2
            ddy1=ddy2; ddy2=ddy;
        end
        if nargout>3
            dddy1=dddy2; dddy2=dddy;
        end
    end
end
end
%% Demo
% DrawBasePyram(5,1,0);
% DrawResGalerkin(18);
%函数体内有参数说明
% 以下附上源码，取消一层注释即可新建函数 DrawBasePyram和DrawResGalerkin
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function DrawBasePyram(m,BaseType,opt)
% % Input:
% % BaseType:可视化基函数的类型（1——顶点基；2——边基；3——面基；4——体基）
% % opt:可视化基函数的导数类型（0——函数值G；1——Gr；2——Gs；3——Gt）
% % m:基函数阶次选择，m建议别选太大，否则基函数数量太多，运行时长太久
% % 体函数由于数量较多，只显示前10个基
% 
% number.edge=[m-1,m-1,m-1,m-1,m-1,m-1,m-1,m-1];
% number.face={[m-1,m-1],m-2,m-2,m-2,m-2};
% number.H=m-2;
% n=50;
% %% 金字塔网格点
% % [r,s,t,tri]=StdPyramPoint(n);
% [R, S, T, ~]=meshgridpyramc(n, n, n);
% r=R(:);s=S(:);t=T(:);
% [G,Gr,Gs,Gt]=BasePyram(r,s,t,number);
% 
% switch BaseType
%     case 1
%         for i=1:5
%             set(0,'defaultfigurecolor','w');
%             figure; hold on;
%             title(['Vertex',num2str(i)]);
%             
%             switch opt
%                 case 0
%                     c=G(:,i);
%                 case 1
%                     c=Gr(:,i);
%                 case 2
%                     c=Gs(:,i);
%                 case 3
%                     c=Gt(:,i);
%             end
%             C=reshape(c, n, n, n);
%             hexasurfg(R, S, T, C, 2);
%             
%             colorbar;
%             view(3); alpha(1);
%             axis equal;
%             axis off
%             shading interp
%             
%             
%         end
%     case 2
%         for i=6:5+8*(m-1)
%             set(0,'defaultfigurecolor','w');
%             figure; hold on;
%             title(['Edge',num2str(i-5)])
%             
%             switch opt
%                 case 0
%                     c=G(:,i);
%                 case 1
%                     c=Gr(:,i);
%                 case 2
%                     c=Gs(:,i);
%                 case 3
%                     c=Gt(:,i);
%             end
%             C=reshape(c, n, n, n);
%             hexasurfg(R, S, T, C, 2);
%             
%             
%             colorbar;
%             view(3);
%             alpha(1);
%             axis equal;
%             axis off
%             shading interp
%         end
%     case 3
%         for i=6+8*(m-1):5+8*(m-1)+(m-1)^2+2*(m-1)*(m-2)
%             set(0,'defaultfigurecolor','w');
%             figure; hold on;
%             title(['Face',num2str(i-6-8*(m-1))])
%             switch opt
%                 case 0
%                     c=G(:,i);
%                 case 1
%                     c=Gr(:,i);
%                 case 2
%                     c=Gs(:,i);
%                 case 3
%                     c=Gt(:,i);
%             end
%             C=reshape(c, n, n, n);
%             hexasurfg(R, S, T, C, 2);
%             colorbar;
%             view(3); alpha(1);
%             axis equal;
%             axis off
%             shading interp
%         end
%     case 4
%         for i=6+8*(m-1)+(m-1)^2+2*(m-1)*(m-2):5+8*(m-1)+(m-1)^2+2*(m-1)*(m-2)+10
%             set(0,'defaultfigurecolor','w');
%             figure; hold on;
%             title(['Body',num2str(i-5-8*(m-1)-(m-1)^2-2*(m-1)*(m-2))])
%             switch opt
%                 case 0
%                     c=G(:,i);
%                 case 1
%                     c=Gr(:,i);
%                 case 2
%                     c=Gs(:,i);
%                 case 3
%                     c=Gt(:,i);
%             end
%             C=reshape(c, n, n, n);
%             hexasurfg(R, S, T, C, 2);
%             colorbar;
%             view(3); alpha(0.1);
%             axis equal;
%             axis off
%             shading interp
%         end
% end
% end
% 
% function [R, S, T, C]=meshgridpyramc(tt, cc, K)
% 
% % Replicates the grid vectors rgv, sgv, tgv and the associated 
% %     integration weights to produce the coordinates of 
% %     a pyramid grid (R, S, T).
% % 
% % Calling Sequences:
% %
% %     [R, S, T, C]=meshgridpyramidc(tt, cc)
% %
% % INPUTS:
% % 
% %      tt   -  If tt is a cell array, tt = {tu, tv, tw} of the parametric coordinates.
% %            All tu, tv and tw should be defined in [0, 1].
% %              If both inputs are nunvers, tt should be the number intergration 
% %              nodes on an edge of (r, s) plane, and cc should be the number 
% %              intergration nodes on t direction.
% %
% %     cc  - integration weights in one dimensional for u, v and w direction 
% %             respectively, cc is a cell {cu, cv, cw}.
% %
% %     K - Int number of nodes on t direction.
% %
% % OUTPUT:
% % 
% %     R    :   a matrix of u coordinates for integration on a pyramid
% %     S    :   a matrix of v coordinates for integration on a pyramid
% %     T    :   a matrix of w coordinates for integration on a pyramid
% %     C    :   a matrix of integration weights for pyramid
% % 
% 
% if nargin==3
%     if isnumeric(tt) && isnumeric(cc)
%         M=tt; N=cc;
%         [r, Cr]=GaussLobattoR(M, -1, 1);
%         [s, Cs]=GaussLobattoR(N, -1, 1);
%         [t, Ct]=GaussLobattoR(K, 0, 1);
%         tt={r, s, t}; cc={Cr, Cs, Ct};
%     end
% end
% 
% [s, r, t]=meshgrid(tt{2}, tt{1}, tt{3}); 
% [Cs, Cr, Ct]=meshgrid(cc{2}, cc{1}, cc{3}); 
% 
% % Get nodes on a unit pyramid
% R=(1-t).*r;
% S=(1-t).*s;
% T=t;
% C=(1-t).^2.*Cr.*Cs.*Ct;
% 
% 
% %% demo
% % % Get nodes on a unit hexahedron
% % n=5;
% % [r, Cr]=GaussLobattoR(n, -1, 1);
% % [s, Cs]=GaussLobattoR(n+1, -1, 1);
% % [t, Ct]=GaussLobattoR(n+2, 0, 1);
% % [R, S, T, C]=meshgridpyramc({r, s, t}, {Cr, Cs, Ct});
% % 
% % % Symbbolic expression
% % syms r s t
% % Eq=1+r^6+s^3+t^2;
% % fh=@(r, s, t) 1+ r.^6+s.^3+t.^2;
% % 
% % iEq=int(Eq, r, t-1, 1-t);
% % iEq=int(iEq, s, t-1, 1-t);
% % iEq=double(int(iEq, t, 0, 1));
% % 
% % % Numerical intergation
% % Z=fh(R, S, T);
% % iZ=sum(C(:).*Z(:));
% % 
% % % Draw the domain
% % hexasurfg(R, S, T);
% end
% 
% function hexasurfg(x, y, z, C, k)
% % 
% % hexasurfg : Plot a hexahedron region by surfaces
% % 
% % Calling Sequences:
% %
% %     hexasurfg(x, y, z)
% %
% %     hexasurfg(x, y, z, C)
% %
% %     hexasurfg(x, y, z, C, k)
% % 
% % INPUTS:
% % 
% %     x, y, z   -    Matrices of parametric coordinates of a unit
% %               hexahedron.  See also  meshgridtetrac.
% %
% %     C - A matrix of color
% %
% %     k  -  Step leangth of plot. For N grid, the surfaces 1:k:N 
% %            of each direction will be plotted.
% %
% 
% if nargin==3
%     C=z; k=1;
% elseif nargin==4
%     k=1;
% end
% 
% hold_flag = ishold;
% hold on;
% 
% M=size(y, 1); N=size(y, 2); K=size(y, 3);
% for i=1:k:M
%     surf(squeeze(x(i,:,:)), squeeze(y(i,:,:)), squeeze(z(i,:,:)), squeeze(C(i,:,:))); 
% end
% for i=1:k:N
%     surf(squeeze(x(:,i,:)), squeeze(y(:,i,:)), squeeze(z(:,i,:)), squeeze(C(:,i,:))); 
% end
% for i=1:k:K
%     surf(squeeze(x(:,:,i)), squeeze(y(:,:,i)), squeeze(z(:,:,i)), squeeze(C(:,:,i))); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.1);
% 
% if (~hold_flag)
%     hold off
% end
% end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function DrawResGalerkin(n)
% % n=18; 
% number.edge=[n-1,n-1,n-1,n-1,n-1,n-1,n-1,n-1];
% number.face={[n-1,n-1],n-2,n-2,n-2,n-2};
% number.H=n-2; M=n+2;
% 
% % Integration nodes
% na=4; nx=M+na; ny=M+na; nz=M+na;
% [x, y, z, c]=meshgridpyramc(nx, ny, nz);
% 
% % Get the basis on the unti tetrahedron
%  [G, Gx, Gy, Gz]=BasePyram(x(:), y(:), z(:), number);
% 
% % Function values and derivatives
% %% Function Test1
% fun=@(x,y,z) cos(x).*cos(y).*cos(z);
% funx=@(x,y,z) -sin(x).*cos(y).*cos(z);
% funy=@(x,y,z) -cos(x).*sin(y).*cos(z);
% funz=@(x,y,z) -cos(x).*cos(y).*sin(z);
% %% Function Test2
% % fun=@(x,y,z) (x.^3).*(y.^3).*(z.^8)+(x.^2).*(y)+2*(x)+(y)+1;
% % funx=@(x,y,z) (3*x.^2).*(y.^3).*(z.^8)+(2*x).*(y)+2;
% % funy=@(x,y,z) (x.^3).*(3*y.^2).*(z.^8)+(x.^2)+1;
% % funz=@(x,y,z) (x.^3).*(y.^3).*(8*z.^7);
% %% Function Test3
% % fun=@(x,y,z) (x.^3)+(z.^8);
% % funx=@(x,y,z) (3*x.^2);
% % funy=@(x,y,z) 0;
% % funz=@(x,y,z) (8*z.^7);
% %% Function Test4
% % fun=@(x,y,z) x.*y./(1.1-z);
% % funx=@(x,y,z) y./(1.1-z);
% % funy=@(x,y,z) x./(1.1-z);
% % funz=@(x,y,z) x.*y./((1.1-z).^2);
% %% Function Test5
% % fun=@(x,y,z) x+y+z.^7+100.*z;
% % funx=@(x,y,z) 0*x+1;
% % funy=@(x,y,z) 0*x+1;
% % funz=@(x,y,z) 7*z.^6+100;
% 
% F=fun(x,y,z); Fx=funx(x,y,z); 
% Fy=funy(x,y,z); Fz=funz(x,y,z); 
% 
% % Galerkin interpolation
% 
% CJ=diag(c(:)); 
% Gk=(G'*CJ*G)\G'*CJ; 
% w=Gk*F(:); 
% 
% % % Plot results
% figure; hexasurfg(x, y, z, F, 2); 
% title('F'); axis off
% figure; hexasurfg(x, y, z, Fx, 2); 
% title('Fx - Exact'); axis off
% figure; hexasurfg(x, y, z, Fy, 2); 
% title('Fy - Exact'); axis off
% figure; hexasurfg(x, y, z, Fz, 2); 
% title('Fz - Exact'); axis off
% 
% Fi=reshape(G*w, nx, ny, nz);
% figure; hexasurfg(x, y, z, F-Fi, 2); 
% title('Res-F-Galerkin');
% axis off
% 
% Fxi=reshape(Gx*w, nx, ny, nz);
% figure; hexasurfg(x, y, z, Fx-Fxi, 2); 
% title('Res-Fx-Galerkin'); 
% axis off
% 
% Fyi=reshape(Gy*w, nx, ny, nz);
% figure; hexasurfg(x, y, z, Fy-Fyi, 2); 
% title('Res-Fy-Galerkin');
% axis off
% 
% Fzi=reshape(Gz*w, nx, ny, nz);
% figure; hexasurfg(x, y, z, Fz-Fzi, 2); 
% title('Res-Fz-Galerkin');
% axis off
% end
% 
% function [R, S, T, C]=meshgridpyramc(tt, cc, K)
% 
% % Replicates the grid vectors rgv, sgv, tgv and the associated 
% %     integration weights to produce the coordinates of 
% %     a pyramid grid (R, S, T).
% % 
% % Calling Sequences:
% %
% %     [R, S, T, C]=meshgridpyramidc(tt, cc)
% %
% % INPUTS:
% % 
% %      tt   -  If tt is a cell array, tt = {tu, tv, tw} of the parametric coordinates.
% %            All tu, tv and tw should be defined in [0, 1].
% %              If both inputs are nunvers, tt should be the number intergration 
% %              nodes on an edge of (r, s) plane, and cc should be the number 
% %              intergration nodes on t direction.
% %
% %     cc  - integration weights in one dimensional for u, v and w direction 
% %             respectively, cc is a cell {cu, cv, cw}.
% %
% %     K - Int number of nodes on t direction.
% %
% % OUTPUT:
% % 
% %     R    :   a matrix of u coordinates for integration on a pyramid
% %     S    :   a matrix of v coordinates for integration on a pyramid
% %     T    :   a matrix of w coordinates for integration on a pyramid
% %     C    :   a matrix of integration weights for pyramid
% % 
% 
% if nargin==3
%     if isnumeric(tt) && isnumeric(cc)
%         M=tt; N=cc;
%         [r, Cr]=GaussLobattoR(M, -1, 1);
%         [s, Cs]=GaussLobattoR(N, -1, 1);
%         [t, Ct]=GaussLobattoR(K, 0, 1);
%         tt={r, s, t}; cc={Cr, Cs, Ct};
%     end
% end
% 
% [s, r, t]=meshgrid(tt{2}, tt{1}, tt{3}); 
% [Cs, Cr, Ct]=meshgrid(cc{2}, cc{1}, cc{3}); 
% 
% % Get nodes on a unit pyramid
% R=(1-t).*r;
% S=(1-t).*s;
% T=t;
% C=(1-t).^2.*Cr.*Cs.*Ct;
% 
% 
% %% demo
% % % Get nodes on a unit hexahedron
% % n=5;
% % [r, Cr]=GaussLobattoR(n, -1, 1);
% % [s, Cs]=GaussLobattoR(n+1, -1, 1);
% % [t, Ct]=GaussLobattoR(n+2, 0, 1);
% % [R, S, T, C]=meshgridpyramc({r, s, t}, {Cr, Cs, Ct});
% % 
% % % Symbbolic expression
% % syms r s t
% % Eq=1+r^6+s^3+t^2;
% % fh=@(r, s, t) 1+ r.^6+s.^3+t.^2;
% % 
% % iEq=int(Eq, r, t-1, 1-t);
% % iEq=int(iEq, s, t-1, 1-t);
% % iEq=double(int(iEq, t, 0, 1));
% % 
% % % Numerical intergation
% % Z=fh(R, S, T);
% % iZ=sum(C(:).*Z(:));
% % 
% % % Draw the domain
% % hexasurfg(R, S, T);
% end
% 
% function hexasurfg(x, y, z, C, k)
% % 
% % hexasurfg : Plot a hexahedron region by surfaces
% % 
% % Calling Sequences:
% %
% %     hexasurfg(x, y, z)
% %
% %     hexasurfg(x, y, z, C)
% %
% %     hexasurfg(x, y, z, C, k)
% % 
% % INPUTS:
% % 
% %     x, y, z   -    Matrices of parametric coordinates of a unit
% %               hexahedron.  See also  meshgridtetrac.
% %
% %     C - A matrix of color
% %
% %     k  -  Step leangth of plot. For N grid, the surfaces 1:k:N 
% %            of each direction will be plotted.
% %
% 
% if nargin==3
%     C=z; k=1;
% elseif nargin==4
%     k=1;
% end
% 
% hold_flag = ishold;
% hold on;
% 
% M=size(y, 1); N=size(y, 2); K=size(y, 3);
% for i=1:k:M
%     surf(squeeze(x(i,:,:)), squeeze(y(i,:,:)), squeeze(z(i,:,:)), squeeze(C(i,:,:))); 
% end
% for i=1:k:N
%     surf(squeeze(x(:,i,:)), squeeze(y(:,i,:)), squeeze(z(:,i,:)), squeeze(C(:,i,:))); 
% end
% for i=1:k:K
%     surf(squeeze(x(:,:,i)), squeeze(y(:,:,i)), squeeze(z(:,:,i)), squeeze(C(:,:,i))); 
% end
% shading interp;
% colorbar;
% view(3); alpha(0.1);
% 
% if (~hold_flag)
%     hold off
% end
% end

















