clear all;

alpha = 0.4;
beta = 0.96;
delta = 0.1;
Abar = 1.0;

kss = (alpha*beta*Abar/(1-beta*(1-delta)))^(1/(1-alpha));
css = kss^alpha - delta*kss;

%T = 30; % > 150 does not work (when delta=0.1)
T = 150; % > 35 does not work (when delta=1.0)

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

save('fig_4_4.mat');

figure;
plot(k(1:100),'b-','LineWidth',2);
hold on;
plot(1,k(1),'ro-','LineWidth',2);
plot(100,k(100),'ro-','LineWidth',2);
% plot([1,100],[kss,kss],'k--');
xlabel('期間','FontSize',14);
ylabel('k','FontSize',14);
%title('k');
xlim([1 100]);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
saveas(gcf,'pf.jpg','jpg')
print -depsc2 pf.eps;

figure;
plot(k,c,'b*','MarkerSize',10);
hold on;
plot(k(1),c(1),'r*','MarkerSize',10);
plot(k(100),c(100),'r*','MarkerSize',10);
xlabel('k','FontSize',14);
ylabel('c','FontSize',14);
box on;
grid on;
set(gca,'Fontsize',14)
% set(gca,'FontName','Times New Roman')
set(gcf,'color','w')
saveas(gcf,'pf2.jpg','jpg')
print -depsc2 pf2.eps;