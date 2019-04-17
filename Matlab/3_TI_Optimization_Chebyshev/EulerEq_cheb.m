function  res = EulerEq_cheb(cons,m,capital,theta)
% cを与えたときのオイラー方程式の残差を返す関数

wealth = capital.^m.alpha + (1.-m.delta).*capital;

kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(m.kgrid(1),kprime);

% 次期の政策関数を線形補間: m.nk=21のときは政策関数の形がおかしい???
%cnext = interp1(m.kgrid,cfcn,kprime,'linear','extrap');
% 次期の価値関数をスプライン補間
%cnext = interp1(m.kgrid,cfcn,kprime,'spline');
% 次期の価値関数を多項式補間
T = polybas(m.kmin,m.kmax,m.nk,kprime);
cnext = T*theta;

% オイラー方程式
res = (1/cons) - m.beta*(1/cnext)*(m.alpha*kprime.^(m.alpha-1) + (1.-m.delta));
 
return