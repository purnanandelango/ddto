function [r,v,p,y,T,s] = disassemble(x,u,n,ntarg_p_1)

    N = size(x,2);

    r = zeros(n,N,ntarg_p_1);
    v = zeros(n,N,ntarg_p_1);
    p = zeros(N,ntarg_p_1);
    y = zeros(N,ntarg_p_1);
    
    T = zeros(n,N,ntarg_p_1);
    s = zeros(N,ntarg_p_1);

    for k = 1:N
    
        xk = reshape(x(:,k),[2*n+2,ntarg_p_1]);
        uk = reshape(u(:,k),[n+1,ntarg_p_1]);
    
        r(:,k,:) = xk(1:n,:);
        v(:,k,:) = xk(n+1:2*n,:);
        p(k,:)   = xk(2*n+1,:);
        y(k,:)   = xk(2*n+2,:);
    
        T(:,k,:) = uk(1:n,:);
        s(k,:)   = uk(n+1,:);
    
    end

end