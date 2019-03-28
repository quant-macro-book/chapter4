function util = CRRA(cons,gamma)
% Function CRRA
%  utility = CRRA( consumption, gamma )
%
% Purpose:
%  Compute CRRA utility function
%
%  Record of revisions:
%     Date     Programmer  Description of change
%  ==========  ==========  =====================
%  10/05/2002  T. Yamada   Original code

if gamma ~= 1;
    util = cons.^(1-gamma)./(1-gamma);
else
    util = log(cons);
end

return;
