function  res = EulerEq_cheb(cons,capital,theta)
% オイラー方程式に代入した際の残差を返す関数

global beta gamma alpha delta kmin kmax nk kgrid

wealth = capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(kgrid(1),kprime);

% 次期の価値関数を多項式補間
% Tはk=kprimeで評価した基底関数
T = polybas(kmin,kmax,nk,kprime);
% 多項式の係数thetaを基底関数に掛けて近似値を求める
cnext = T*theta;

% オイラー方程式の残差を求める（u'(c)をmu_CRRA関数を用いて計算している）
%res = mu_CRRA(cons,gamma) - beta*mu_CRRA(cnext,gamma)*(alpha*kprime.^(alpha-1) + (1.-delta));
res = (1/cons) - beta*(1/cnext)*(alpha*kprime.^(alpha-1) + (1.-delta));
 
return