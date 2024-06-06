function err = calcerr(invT,cfcn0)
%% �I�C���[����������덷�𑪒�

global beta gamma alpha delta kmin kmax nk kgrid

% ���̃O���b�h�ł̓I�C���[�������̌덷�̓[���ɂȂ邽�߁A�O���b�h���ׂ����Ƃ�
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