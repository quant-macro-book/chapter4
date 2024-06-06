clear all;

alpha = 0.33; %1/3;
beta = 0.95;
delta = 0.5; %1.0; %0.1;
gamma = 2.0;
Abar = 1.0;
gbar = 0.2;

% % LS
% alpha = 0.33; %1/3;
% beta = 0.95;
% delta = 0.2; %0.1;
% gamma = 2.0;
% Abar = 1.0;
% gbar = 0.2;

% # steady state (all taus are zero)
% kss = (α*β*Abar/(1-β*(1-δ)))**(1/(1-α))
% css = Abar*kss**α - δ*kss - gbar
% Rss = 1-δ+α*Abar*kss**(α-1)
% print([kss,css])
kss = (alpha*beta*Abar/(1-beta*(1-delta)))^(1/(1-alpha));
css = Abar*kss^alpha - delta*kss - gbar;
Rss = 1-delta+alpha*Abar*kss^(alpha-1);

T = 1000;

%k0 = 0.1*kss;
k0 = kss;

% imaginary state
kvec0 = linspace(k0,kss,T)';
kvec1 = kvec0;

tauc = zeros(T,1);
%tauc(10:T) = 0.2; 
tauc(10:T) = 0.2; 

cvec = zeros(T,1);
kvec = zeros(T,1);
rvec = zeros(T,1);

damp = 0.9;
diff = 1e+4;
crit = 1e-4;
iter = 0;

while(diff>crit)

    rvec = alpha*Abar*kvec0.^(alpha-1);
    wvec = (1-alpha)*Abar*kvec0.^alpha;
        
    cvec(T) = css;
    for i=T-1:-1:1 % t = T-1,...1,0
        cvec(i) = (1/beta/(1+rvec(i+1)-delta)*(1+tauc(i+1))/(1+tauc(i))*(cvec(i+1)^gamma))^(1/gamma);
    end
    
    kvec(1) = k0;
    for j=1:1:T-1
%        kvec(j+1) = Abar*kvec0(j)^alpha + (1-delta)*kvec(j) - gbar - cvec(j); % k1,k2,...,kT
        kvec(j+1) = rvec(j)*kvec0(j) + wvec(j) + (1-delta)*kvec(j) - gbar - cvec(j); % k1,k2,...,kT
    end
    
    % update
    kvec1 = damp*kvec0 + (1-damp)*kvec;
    
    diff = max(abs(kvec1-kvec0));
    kvec0 = kvec1;
    iter = iter + 1;

    % use the terminal condition as the criterion
    disp([iter diff]);
        
end

figure;
subplot(241);
plot(kvec);
xlim([1 40]);
subplot(242);
plot(cvec);
xlim([1 40]);
subplot(243);
plot(rvec+1-delta);
xlim([1 40]);
subplot(244);
plot(tauc,'r-');
xlim([1 40]);