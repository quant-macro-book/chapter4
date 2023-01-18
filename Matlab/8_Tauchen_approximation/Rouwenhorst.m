function [Z, PI] = Rouwenhorst(rho, sigmas, znum)

% Creates a discrete approximation to a first order autoregressive process
% with serial correlation p + q - 1.  See Rouwenhorst (1995) (pps. 325-329)
% in Cooley <i>Frontiers of Business Cycle Research</i> Princeton.

if (nargin < 3)	% default is 2x2 transition matrix
	znum = 2;
end 

p = (rho + 1)/2;
q = p;

hlag = 1; h = [];

for i = 2:1:znum
	vec0 = zeros(i-1,1);
	h = p*[hlag vec0; vec0' 0] + (1-p)*[vec0 hlag; 0 vec0'] + ...
	 	(1- q)*[vec0' 0; hlag vec0] + q*[0 vec0'; vec0 hlag];
	h(2:i-1,:) = h(2:i-1,:)./2; hlag = h;
	
end

PI = h;

% symmetrically and evenly spaced between [-epsilon, epsilon] with h elements.  
% When p = q, then then variance of shock is epsilon^2/(h-1).  
zvar = (sigmas^2)/( 1 - rho^2);
epsilon = sqrt((znum - 1)*zvar);

Z = linspace(-epsilon, epsilon, znum)';