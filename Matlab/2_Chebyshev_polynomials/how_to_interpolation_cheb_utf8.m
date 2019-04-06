%% 
close all;

%% DATA POINTS
xmin = -1;
xmax = 1;
nxd = 21;
xd = linspace(xmin,xmax,nxd)';
yd = 1./(1+25*xd.^2);
nx = 1001;
x0 = linspace(xmin,xmax,nx)';
y0 = 1./(1+25*x0.^2);

figure;
plot(x0, y0, 'k-');
hold on;
plot(xd, yd, 'o', 'color', 'blue', 'MarkerSize', 12, 'linewidth', 3);
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig_data.eps','epsc2');

%% MATLAB関数を使った内挿

x1 = linspace(xmin,xmax,nx)';
y1 = zeros(nx, 1);

for i = 1:nx
    y1(i, 1) = interp1(xd, yd, x1(i), 'linear', 'extrap');
end

%% ordinary polynomial
X = ones(nxd,1);
X2 = ones(nx,1);
x2 = x1;
for i = 1:nxd-1
    X = [X xd.^i];
    X2 = [X2 x2.^i];
end
% ordinary least squares
b = (X'*X)\(X'*yd);
y2 = X2*b;

figure;
plot(x1, y1, '-', 'color', 'blue', 'linewidth', 3);
hold on;
plot(x2, y2, '--', 'color', 'red', 'linewidth', 3);
plot(x0, y0, 'k-');
legend('線形近似', '多項式近似', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig_interp.eps','epsc2');

%% Chebyshev polynomial
xcheb = polygrid(xmin,xmax,nxd);
ycheb = 1./(1+25*xcheb.^2);
T = polybas(xmin,xmax,nxd,xcheb);
theta = T\ycheb;

x3 = x1;
T3 = polybas(xmin,xmax,nxd,x3);
% fitting polynomial
y3 = T3*theta;

figure;
plot(x3, y3, '-', 'color', 'blue', 'linewidth', 3); hold('on');
hold on;
plot(xcheb, ycheb, '*', 'color', 'blue', 'MarkerSize',12, 'linewidth', 3); hold('on');
plot(x0, y0, 'k-');
legend('多項式近似','評価点','Location', 'NorthEast');
grid on;
set(gca,'Fontsize',16);
saveas (gcf,'Fig_cheb_n21.eps','epsc2');

return;
