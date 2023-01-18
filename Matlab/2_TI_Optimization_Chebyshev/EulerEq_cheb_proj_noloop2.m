function res = EulerEq_cheb_proj_noloop2(theta_new,T0,theta)

global beta gamma alpha delta kmin kmax nk kgrid

res_vec = zeros(nk,1);

capital = kgrid;
wealth = capital.^alpha + (1.-delta).*capital;
cons = T0*theta_new;
kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(kgrid(1),kprime);        

% 次期の価値関数を多項式補間
% Tはk=kprimeで評価した基底関数
T1 = polybas(kmin,kmax,nk,kprime);
% 多項式の係数thetaを基底関数に掛けて近似値を求める
if (nargin == 2)
    cnext = T1*theta_new; % 新しいthetaの値を用いる
else
    cnext = T1*theta; % 古いthetaの値を用いる
end

res_vec = (1./cons) - beta*(1./cnext).*(alpha*kprime.^(alpha-1) + (1.-delta));
res = sum(res_vec.^2);

end