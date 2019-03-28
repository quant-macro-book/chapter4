function T = polybas(kmin,kmax,Np,kgrid)

% Np: the order of polynomial+1, Ng: # of grid points
Ng = length(kgrid);
% transform from k to x
x = 2.*(kgrid-kmin)/(kmax-kmin) - ones(Ng,1);

T = zeros(Ng,Np);
T0 = ones(Ng,1);
T1 = x;
T2 = 2.*x.*T1 - T0;
T(:,1) = T1;
T(:,2) = T2;

for i=3:Np
    T(:,i) = 2.*x.*T(:,i-1) - T(:,i-2);
end

T = [T0 T(:,1:(Np-1))];