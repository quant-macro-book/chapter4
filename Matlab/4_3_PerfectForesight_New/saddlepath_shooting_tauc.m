clear all;

alpha = 0.33; %1/3;
beta = 0.95;
delta = 0.2;
gamma = 2.0
Abar = 1.0;
gbar = 0.2;

kss = (alpha*beta*Abar/(1-beta*(1-delta)))^(1/(1-alpha));
css = kss^alpha - delta*kss - gbar;

T = 100;

k0 = kss;

% % 1. using fzero
% c0 = 0.129139420765888;
% options = optimset('TolX',1e-16);
% c0 = fzero(@shooting_err,c0,options,k0,T,alpha,beta,delta,Abar,kss);
% [k c] = shooting(c0,k0,T,alpha,beta,delta,Abar);

tauc = zeros(T+1,1);
tauc(11:T) = 0.2;

% 2. bisection
cL = 0.001;
cH = 3.0;
% for i = 1:100
diff = 1e+4;
while(diff>1e-15)

    c0 = (cL+cH)/2;
    [k c R] = shooting_tauc(c0,k0,T,alpha,beta,delta,gamma,Abar,gbar,tauc);
    err = k(T)-kss;
    disp([c0 err]);
    if err>0
        cL = c0;
    else
        cH = c0;
    end

    diff = abs(cH-cL);

end

figure;
subplot(241);
plot(k,'b-','LineWidth',2.0);
title('k');
xlim([1 40]);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(242);
plot(c,'b-','LineWidth',2.0);
title('c');
xlim([1 40]);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(243);
plot(R,'b-','LineWidth',2.0);
title('R');
xlim([1 40]);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
subplot(244);
plot(tauc,'r-','LineWidth',2.0);
title('\tau_{c}');
xlim([1 40]);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
saveas(gcf,'PF_tauc.jpg','jpg')
print -depsc2 PF_tauc.eps
%plot(k);
%xlim([1 40]);
%plot(k,c,'*');