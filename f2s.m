function [ ii_ ] = f2s( L, nd)


    ii_=zeros(nd^3,1);
    
    for n=1:size(L,1)
%         i = ii(n);
        Li = L(n,:);
        [l1_,maxi] = max(Li);
        Li(maxi) = [];
        l2_ = Li(1);
        l3_ = Li(2);
        n_ = l1_ + nd*(l2_-1) + nd^2*(l3_-1);
        ii_(round(n_))=1;
    end
    
    

end

