function res = EulerEq(cons,capital,cfcn)
% オイラー方程式に代入した際の残差を返す関数

global beta gamma alpha delta kgrid

wealth = capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(kgrid(1),kprime);

% 次期の政策関数を線形補間
cnext = interp1(kgrid,cfcn,kprime,'linear','extrap');
% 次期の価値関数をスプライン補間
%cnext = interp1(kgrid,cfcn,kprime,'spline');

% オイラー方程式の残差を求める（u'(c)をmu_CRRA関数を用いて計算している）
res = mu_CRRA(cons,gamma) - beta*mu_CRRA(cnext,gamma)*(alpha*kprime.^(alpha-1) + (1.-delta));
 
return