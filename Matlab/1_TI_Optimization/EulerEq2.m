function  res = EulerEq2(cons,m,capital,cfcn)
% Function EulerEq
%  [res] = EulerEq(kprime)
%
% 目的:
% cを与えたときのオイラー方程式の残差を返す関数.
% main_ti.mから呼び出して使う.

wealth = capital.^m.alpha + (1.-m.delta).*capital;

kprime = wealth - cons;
%% トリック: k'は正の値しか取らない
kprime = max(m.kgrid(1),kprime);

% 次期の政策関数を線形補間: m.nk=21のときは政策関数の形がおかしい
%cnext = interp1(m.kgrid,cfcn,kprime,'linear','extrap');

% 次期の価値関数をスプライン補間
cnext = interp1(m.kgrid,cfcn,kprime,'spline');

%% オイラー方程式
res = (1/cons) - m.beta*(1/cnext)*(m.alpha*kprime.^(m.alpha-1) + (1.-m.delta));
 
return