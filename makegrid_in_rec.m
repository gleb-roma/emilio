function [d_new,b_new] = makegrid_in_rec(c,R,K)


    mind=c(1)-R(1); maxd=c(1)+R(1);
    minb=c(2)-R(2); maxb=c(2)+R(2);
    dstep=(maxd-mind)/(K+1);
    bstep=(maxb-minb)/(K+1);
    num_new_=0;
    for i=1:K
        for j=1:K
            d_ = mind + dstep*i;
            b_ = minb + bstep*j;
%             if sum(([d_ b_]-c).^2./R.^2) <= 1
                num_new_=num_new_+1;
                d_new(num_new_,1) = d_;
                b_new(num_new_,1) = b_;
%             end
        end
    end
end

