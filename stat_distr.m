function [ gp_stat, error_code ] = stat_distr( gp )
    error_code=0;
    if norm(sum(gp,2)-ones(size(gp,1),1))>1e-4, warning('gp is not a prob matrix'), error_code=3; return; end
    opts.isreal = 1;
    n=size(gp,1);
    b=gp';
    opts.v0=ones(n,1)/n;
    opts.maxit=500;
    eps=1e-4;
    opts.tol=eps/n;
    [eigve,eigva_]=eigs(b,1,'lm',opts);
    eigva=diag(eigva_);
    [eigva,ee]=sort(eigva,'descend');
    if abs(eigva(1)-1)>eps, warning('Largest eigenvalue neq 1'), error_code=1, end;
    ev = eigve(:,ee(1));
    gp_stat = ev./sum(ev);

    if any(gp_stat<-eps/n), warning('gp_stat is not a positive measure'), end
    % delete negative entries and re normalize
    ii=(gp_stat<0);
    gp_stat(ii)=0;  
    if sum(gp_stat)>1.01, warning('gp_stat is TOO non-positive'), error_code=2; end;
    gp_stat = gp_stat./sum(gp_stat);


end

