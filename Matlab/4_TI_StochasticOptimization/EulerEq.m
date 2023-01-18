function res = EulerEq(cons,capital,cfcn)
% �I�C���[�������ɑ�������ۂ̎c����Ԃ��֐�

global beta gamma alpha delta kgrid

wealth = capital.^alpha + (1.-delta).*capital;

kprime = wealth - cons;
% �g���b�N: k'�͐��̒l�������Ȃ�
kprime = max(kgrid(1),kprime);

% �����̐���֐�����`���
cnext = interp1(kgrid,cfcn,kprime,'linear','extrap');
% �����̉��l�֐����X�v���C�����
%cnext = interp1(kgrid,cfcn,kprime,'spline');

% �I�C���[�������̎c�������߂�iu'(c)��mu_CRRA�֐���p���Čv�Z���Ă���j
res = mu_CRRA(cons,gamma) - beta*mu_CRRA(cnext,gamma)*(alpha*kprime.^(alpha-1) + (1.-delta));
 
return