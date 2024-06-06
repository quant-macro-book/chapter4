function res_vec = EulerEq_cheb_proj(theta_new,T0,theta)

global beta gamma alpha delta kmin kmax nk kgrid

res_vec = zeros(nk,1);

% ループは必要ない？
for i = 1:nk

    capital = kgrid(i);
    wealth = capital.^alpha + (1.-delta).*capital;
%    T0 = polybas(kmin,kmax,nk,capital);
%    cons = T0*theta_new;
    cons = T0(i,:)*theta_new;

    kprime = wealth - cons;
    % トリック: k'は正の値しか取らない
    kprime = max(kgrid(1),kprime);        

    % 次期の価値関数を多項式補間
    % Tはk=kprimeで評価した基底関数
    T1 = polybas(kmin,kmax,nk,kprime);
    % 多項式の係数thetaを基底関数に掛けて近似値を求める
%    cnext = T1*theta; % 古いthetaの値を用いる
    cnext = T1*theta_new; % 新しいthetaの値を用いる

    % オイラー方程式の残差を求める（u'(c)をmu_CRRA関数を用いて計算している）
    %res = mu_CRRA(cons,gamma) - beta*mu_CRRA(cnext,gamma)*(alpha*kprime.^(alpha-1) + (1.-delta));
    res_vec(i) = (1/cons) - beta*(1/cnext)*(alpha*kprime.^(alpha-1) + (1.-delta));

end

end