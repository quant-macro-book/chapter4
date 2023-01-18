function [Z,Zprob] = tauchen(N,mu,rho,sigma,m)

Z     = zeros(N,1); % グリッド
Zprob = zeros(N,N); % 遷移確率の行列
c     = (1-rho)*mu; % 定数項

% 等間隔のグリッドを定める
% 最大値と最小値
zmax  = m*sqrt(sigma^2/(1-rho^2));
zmin  = -zmax;
% グリッド間の間隔
w = (zmax-zmin)/(N-1);

Z = linspace(zmin,zmax,N)';
% 定常状態はmu
Z = Z + mu;

% グリッドを所与として、遷移確率を求める
for j = 1:N
    for k = 1:N
        if k == 1
            Zprob(j,k) = cdf_normal((Z(1)-c-rho*Z(j)+w/2)/sigma);
        elseif k == N
            Zprob(j,k) = 1 - cdf_normal((Z(N)-c-rho*Z(j)-w/2)/sigma);
        else
            Zprob(j,k) = cdf_normal((Z(k)-c-rho*Z(j)+w/2)/sigma) - ...
                         cdf_normal((Z(k)-c-rho*Z(j)-w/2)/sigma);
        end
    end
end

% 正規分布の累積分布関数
function c = cdf_normal(x)
    c = 0.5*erfc(-x/sqrt(2));