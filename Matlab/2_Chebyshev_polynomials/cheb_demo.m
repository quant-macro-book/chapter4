clear all;

x = linspace(-1,1,101)';

% ordinary polynomial
figure;
plot(x,x.^2,'k-','LineWidth',1.0);
hold on;
plot(x,x.^4,'k--','LineWidth',1.0);
plot(x,x.^8,'k:','LineWidth',1.0);
legend('x^2','x^4','x^8');
axis square;
set(gca,'Fontsize',16);
saveas (gcf,'ordipoly.eps','epsc2');
%print -depsc2 ordipoly.eps

% Chebyshev polynomial
T = polybas(-1,1,9,x);

figure;
plot(x,T(:,3),'k-','LineWidth',1.0);
hold on;
plot(x,T(:,5),'k--','LineWidth',1.0);
plot(x,T(:,9),'k:','LineWidth',1.0);
plot([-1 1],[0 0],'r-');
legend('T_2(x)','T_4(x)','T_8(x)');
axis square;
set(gca,'Fontsize',16);
saveas (gcf,'chebpoly.eps','epsc2');
%print -depsc2 chebpoly.eps