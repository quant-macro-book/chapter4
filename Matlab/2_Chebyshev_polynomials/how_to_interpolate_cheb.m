% 3節：多項式近似
close all;

%% データ点
xmin = -1;
xmax = 1;
nxd = 11;
xd = linspace(xmin, xmax, nxd)'; %collect(LinRange(xmin, xmax, nxd))
yd = f(xd);

%% 関数による値
nx = 1001;
x0 = linspace(xmin, xmax, nx)'; %collect(LinRange(xmin, xmax, nx))
y0 = f(x0);

%%
figure;
plot(x0,y0,'k-','LineWidth',1.5);
hold on;
plot(xd,yd,'bo','MarkerSize',12,'LineWidth',2.0);
grid on;
legend({'$\frac{1}{1+x^2}$'},'Interpreter','latex');
set(gca,'Fontsize',16);
saveas (gcf,'Fig_data.eps','epsc2');

%% Matlab関数(interp1)を使った線形補間による近似
x1 = linspace(xmin, xmax, nx)';
y1 = interp1(xd,yd,x1,'linear','extrap');

%% 通常の多項式による近似
Xd = ones(nxd,nxd);
X2 = ones(nx,nxd);
x2 = x1;
for i = 1:nxd-1
    Xd(:,i+1) = xd.^i;
    X2(:,i+1) = x2.^i;
end

b = (Xd'*Xd)\(Xd'*yd);
y2 = X2*b;

%%
figure;
plot(x0,y0,'k-','LineWidth',1.5);
hold on;
plot(x1,y1,'-','Color','blue','LineWidth',2.0);
plot(x2,y2,'--','Color','red','LineWidth',2.0);
legend({'$\frac{1}{1+x^2}$','線形近似','多項式近似'},'Interpreter','latex');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig_interp.eps','epsc2');

%% チェビシェフ多項式による近似
% N=11
nxd = 11;
xcheb = polygrid(xmin,xmax,nxd);
ycheb = f(xcheb); %#ones(nxd)./(ones(nxd)+25*xcheb.^2)
T = polybas(xmin,xmax,nxd,xcheb);
theta = T\ycheb;

x3 = x1;
T3 = polybas(xmin,xmax,nxd,x3);
y3 = T3*theta;

figure;
plot(x0,y0,'k-','LineWidth',1.5);
hold on;
plot(x3,y3,'-','Color','blue','LineWidth',2.0);
plot(xcheb,ycheb,'*','Color','blue','MarkerSize',12,'LineWidth',2.0);
legend({'$\frac{1}{1+x^2}$','多項式近似','評価点'},'Interpreter','latex'); %,'Location','NorthEast');
grid on;
set(gca,'FontSize',16);
saveas (gcf,'Fig_cheb_n11.eps','epsc2');

% N=21
nxd = 21;
xcheb = polygrid(xmin,xmax,nxd);
ycheb = f(xcheb); %#ones(nxd)./(ones(nxd)+25*xcheb.^2)
T = polybas(xmin,xmax,nxd,xcheb);
theta = T\ycheb;

x3 = x1;
T3 = polybas(xmin,xmax,nxd,x3);
y3 = T3*theta;

figure;
plot(x0,y0,'k-','LineWidth',1.5);
hold on;
plot(x3,y3,'-','Color','blue','LineWidth',2.0);
plot(xcheb,ycheb,'*','Color','blue','MarkerSize',12,'LineWidth',2.0);
legend({'$\frac{1}{1+x^2}$','多項式近似','評価点'},'Interpreter','latex'); %,'Location','NorthEast');
grid on;
set(gca,'FontSize',16);
saveas (gcf,'Fig_cheb_n21.eps','epsc2');