clear all;

alpha = 0.4; %1/3;
beta = 0.96;
delta = 1.0; %0.1; %0.1;
Abar = 1.0;

kss = (alpha*beta*Abar/(1-beta*(1-delta)))^(1/(1-alpha));
css = kss^alpha - delta*kss;

T = 30; % > 150 does not work (when delta=0.1)
%T = 150; % > 35 does not work (when delta=1.0)

k0 = 0.1; %*kss;
% c0 = 0.304329368748233; % solved by Dynare
% [k c] = shooting(c0,k0,T,alpha,beta,delta,Abar);

% 1. using fzero
%c0 = 0.1; %0.304329364503311; %0.304329368748233; % found by Dynare
c0 = 0.129139420765888;
options = optimset('TolX',1e-16);
c0 = fzero(@shooting_err,c0,options,k0,T,alpha,beta,delta,Abar,kss);
[k c] = shooting(c0,k0,T,alpha,beta,delta,Abar);

% % 2. bisection
% cL = 0.29;
% cH = 0.31;
% % for i = 1:100
% diff = 1e+4;
% while(diff>1e-15)
% 
%     c0 = (cL+cH)/2;
%     [k c] = shooting(c0,k0,T,alpha,beta,delta,Abar);
%     err = k(T)-kss;
%     disp([c0 err]);
%     if err>0
%         cL = c0;
%     else
%         cH = c0;
%     end
% 
%     diff = abs(cH-cL);
% 
% end

% k = zeros(T+1,1);
% c = zeros(T+1,1);
% 
% k(1) = k0;
% c(1) = 0.31; %0.304329368748233;
% 
% for t = 1:T
% 
%     k(t+1) = Abar*k(t)^alpha + (1-delta)*k(t) - c(t);
%     c(t+1) = beta*c(t)*(1+Abar*alpha*k(t+1)^(alpha-1)-delta);
% 
% end

figure;
plot(k(1:100));
hold on;
plot(1,k(1),'r*');
plot([1,100],[kss,kss],'k--');
xlabel('time');
ylabel('k');
%title('k');
xlim([1 100]);

%print -depsc2 pf.eps;

figure;
plot(k,c,'*');
hold on;
plot(k(1),c(1),'r*');
xlabel('k');
ylabel('c');

%print -depsc2 pf2.eps;
