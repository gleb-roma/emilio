function [d_new,b_new] = makegrid_in_round(c,R,K)
%     K=5;
%     d_new=zeros(K^2,1); b_new=zeros(K^2,1); 

%     d1=v1(1); d2=v2(1);
%     b1=v1(2); b2=v2(2);
%         R = norm(v1-v2,2);

    mind=c(1)-R(1); maxd=c(1)+R(1);
    minb=c(2)-R(2); maxb=c(2)+R(2);
%     c = (v1+v2)/2;
    dstep=(maxd-mind)/(K+1);
    bstep=(maxb-minb)/(K+1);
    num_new_=0;
    for i=1:K
        for j=1:K
            d_ = mind + dstep*i;
            b_ = minb + bstep*j;
            if sum(([d_ b_]-c).^2./R.^2) <= 1
                num_new_=num_new_+1;
                d_new(num_new_,1) = d_;
                b_new(num_new_,1) = b_;
            end
        end
    end
end

