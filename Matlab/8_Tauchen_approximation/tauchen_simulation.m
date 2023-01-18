% Markov chain simulation

N = 11;
mu = 0;
rho = 0.95;
sigma = 0.01;
m = 2.575; % covers 99% of distribution

[Z,Pz] = tauchen(N,mu,rho,sigma,m)
%[Z,Pz] = Rouwenhorst(rho,sigma,N)

% simulation
% generate exogenous shock sequence
T = 1000;
ivec = zeros(T+1,1);
ivec(1) = ceil(N/2);
% 条件付き累積密度関数：i列の要素の累積和
cumP = cumsum(Pz')';

rng(0);

for t = 1:T
    cumPi = cumP(ivec(t),:);
    % 一様分布から[0,1]の値をとる乱数を発生させ、条件付き累積密度関数と比較する
    % たとえば、rand<cumPi(1)のとき、ivec(t+1)=1となる
    ivec(t+1) = sum(rand-cumPi >= 0)+1;
    ivec(t+1) = min(ivec(t+1),N);
end

figure
plot([0:T],Z(ivec));

% recovering the original AR(1) process
% X = Z(ivec(1:end-1));
% y = Z(ivec(2:end));
% Beta = (X'*X)\X'*y
% e = y-X*Beta;
% s = e'*e/(N-1);
% Sigma = sqrt(s)
% SEBeta = Sigma*sqrt(inv(X'*X))

