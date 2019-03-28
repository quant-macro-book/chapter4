function mu = mu_CRRA( cons, gamma )
% Function mu_CRRA
%  marginal_utility = mu_CRRA( consumption, gamma )
%
% Purpose:
%  Compute marginal utility of CRRA-type function
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  02/22/2016  T. Yamada   Original code

mu = cons.^-gamma;

return;
