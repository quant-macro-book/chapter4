function res = EulerEq(cons,wealth,cfcn)
% オイラー方程式に代入した際の残差を返す関数

global beta gamma alpha delta kgrid zgrid Pzvec

%wealth = technology*capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
% トリック: k'は正の値しか取らない
kprime = max(kgrid(1),kprime);

% 右辺の期待値を計算
% A_{t+1}=A_{j} for j=1,...,nzのそれぞれについて、条件付き確率で加重平均をとることで期待値を計算する
RHS = 0;
for jz = 1:size(Pzvec,1)
    
    % 次期の政策関数を線形補間
    cnext = interp1(kgrid,cfcn(:,jz),kprime,'linear','extrap');
    % 次期の価値関数をスプライン補間
    %cnext = interp1(kgrid,cfcn,kprime,'spline');
    % 条件付き確率で加重平均
    RHS = RHS + Pzvec(jz)*mu_CRRA(cnext,gamma)*(alpha*exp(zgrid(jz))*kprime.^(alpha-1) + (1.-delta));
    
end

% オイラー方程式の残差を求める（u'(c)をmu_CRRA関数を用いて計算している）
res = mu_CRRA(cons,gamma) - beta*RHS;
 
return