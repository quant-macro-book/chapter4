function err = calcerr(invT,cfcn0)
%% オイラー方程式から誤差を測定

global beta gamma alpha delta kmin kmax nk kgrid

% 元のグリッドではオイラー方程式の誤差はゼロになるため、グリッドを細かくとる
theta = invT*cfcn0;
kgrid_err = linspace(kmin,kmax,(nk-1)*10+1)';
T = polybas(kmin,kmax,nk,kgrid_err);
cons = T*theta;
LHS  = mu_CRRA(cons, gamma);

kp   = kgrid_err.^alpha + (1-delta)*kgrid_err - cons;
T = polybas(kmin,kmax,nk,kp);
cnext = T*theta;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cnext,gamma);

err  = RHS./LHS-1.0;