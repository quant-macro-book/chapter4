function [k,c] = shooting(c0,k0,T,alpha,beta,delta,Abar)

k = zeros(T+1,1);
c = zeros(T+1,1);

k(1) = k0;
c(1) = c0;

for t = 1:T

    k(t+1) = max(Abar*k(t)^alpha + (1-delta)*k(t) - c(t),1e-4); % k > 0
    c(t+1) = beta*c(t)*(1+Abar*alpha*k(t+1)^(alpha-1)-delta);

end