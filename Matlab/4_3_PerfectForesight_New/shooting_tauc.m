function [k,c,R] = shooting_tauc(c0,k0,T,alpha,beta,delta,gamma,Abar,gbar,tauc)

k = zeros(T+1,1);
c = zeros(T+1,1);
R = zeros(T+1,1);

k(1) = k0;
c(1) = c0;
R(1) = 1/beta;

for t = 1:T

    k(t+1) = max(Abar*k(t)^alpha + (1-delta)*k(t) - gbar - c(t),1e-4); % k > 0
    R(t+1) = ((1+tauc(t))/(1+tauc(t+1))) * ( 1-delta + alpha*Abar*k(t+1)^(alpha-1) );
    RHS = beta*R(t+1)*(c(t)^gamma);
%    RHS = beta*((1+tauc(t))/(1+tauc(t+1))) * ( 1-delta + alpha*Abar*k(t+1)^(alpha-1) )*(c(t)^gamma);
    c(t+1) = RHS^(1/gamma);

end