function [ ii_ ] = opp_shirks( L, nd, on)

% on is opponent index, 2 or 3 only.
% returns the state indices where %on shirks if %on is freshman
% ii_ is the logical array

    ii_=zeros(nd^3,1);
    
    if on==2
        for n=1:size(L,1)
            Li = L(n,:);
            l1_ = Li(2);
            l2_ = Li(1);
            l3_ = Li(3);
            n_ = l1_ + nd*(l2_-1) + nd^2*(l3_-1);
            ii_(round(n_))=1;

            l1_ = Li(3);
            l2_ = Li(1);
            l3_ = Li(2);
            n_ = l1_ + nd*(l2_-1) + nd^2*(l3_-1);
            ii_(round(n_))=1;
        end
    elseif on==3
        for n=1:size(L,1)
            Li = L(n,:);
            l1_ = Li(2);
            l2_ = Li(3);
            l3_ = Li(1);
            n_ = l1_ + nd*(l2_-1) + nd^2*(l3_-1);
            ii_(round(n_))=1;

            l1_ = Li(3);
            l2_ = Li(2);
            l3_ = Li(1);
            n_ = l1_ + nd*(l2_-1) + nd^2*(l3_-1);
            ii_(round(n_))=1;
        end
    else
        error('on is 2 or 3 only');
    end
    
    

end

