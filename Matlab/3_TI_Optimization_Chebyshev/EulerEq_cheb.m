function  res = EulerEq(cons)
% Function EulerEq
%  [res] = EulerEq(kprime)
%
% 目的:
% k'を一つ与えたときのオイラー方程式の残差を返す関数.
% main_ti.mから呼び出して使う.
%
% グローバル変数: beta gamma alpha delta capital vfcn kgrid

global beta gamma alpha delta capital kmin kmax nk theta %cfcn kgrid

%% オイラー方程式

wealth = capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
%% トリック(1): k'は正の値しか取らない
kprime = max(kmin,kprime);

% 次期の政策関数を線形補間
%cnext = interp1(kgrid, cfcn, kprime, 'linear', 'extrap');

% 次期の価値関数をスプライン補間
%cnext = interp1(kgrid, cfcn, kprime, 'spline');

% 次期の価値関数を多項式補間(Chebyshev)
T = polybas(kmin,kmax,nk,kprime);
cnext = T*theta;

res = (1/cons) - beta*(1/cnext)*(alpha*kprime.^(alpha-1) + (1.-delta));
 
return
