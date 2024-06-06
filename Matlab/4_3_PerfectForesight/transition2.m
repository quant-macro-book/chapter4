clear all;

alpha = 0.4; %1/3;
beta = 0.96;
delta = 1.0; %0.1;
Abar = 1.0;

kss = (alpha*beta*Abar/(1-beta*(1-delta)))^(1/(1-alpha));
css = Abar*kss^alpha - delta*kss;

T = 100;

k0 = 0.1*kss;

% imaginary state
kvec0 = linspace(k0,kss,T)';
kvec1 = kvec0;

yvec = zeros(T,1);
rvec = zeros(T,1);
wvec = zeros(T,1);
cvec = zeros(T,1);
kvec = zeros(T,1);
%evec = zeros(T,1);

lamk = 0.9;
%lame = 0.9;
diff = 1e+4;
crit = 1e-4;
iter = 0;
m    = figure;

while(diff>crit)
% for i = 1:250
    
    yvec = Abar*kvec0.^alpha;
    rvec = alpha*yvec./kvec0;
    wvec = (1-alpha)*yvec;
        
    cvec(T) = css;
    for i=T-1:-1:1 % t = T-1,...1,0
        cvec(i) = 1/beta/(1+rvec(i+1)-delta)*cvec(i+1);
    end
    
    kvec(1) = k0;
    for j=1:1:T-1
        kvec(j+1) = yvec(j) - cvec(j) + (1-delta)*kvec(j); % k1,k2,...,kT
%        kvec(j+1) = Abar*kvec(j)^alpha - cvec(j) + (1-delta)*kvec(j); % k1,k2,...,kT
    end
    
%    evec = ((1-THETA)/ALPHA)*yvec./cvec;
%     evec = (1./h).*(((1-THETA)/ALPHA).*(h./cvec)).^(1/THETA).*kvec0;
    
    % update
    kvec1 = lamk*kvec0 + (1-lamk)*kvec;
%    evec1 = lame*evec0 + (1-lame)*evec;
    
    diff = max(abs(kvec1-kvec0));
%    ediff = max(abs(evec1-evec0));
    kvec0 = kvec1;
%    evec0 = evec1;
    iter = iter + 1;

    % use the terminal condition as the criterion
    disp([iter diff]);
    
%     figure(m);
%     subplot(211);
%     plot([1:T]',kvec0);
%     subplot(212);
%     plot([1:T]',evec0);
    
end

figure;
plot(kvec0);
xlim([1 100]);
