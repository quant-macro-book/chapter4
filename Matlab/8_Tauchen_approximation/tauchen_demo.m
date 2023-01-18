clear all;

N = 11;
mu = 0;
rho = 0.8;
sigma = 1.0;
m = 3.0;

c     = (1-rho)*mu; % 定数項

% 等間隔のグリッドを定める
% 最大値と最小値
zmax  = m*sqrt(sigma^2/(1-rho^2));
zmin  = -zmax;
% グリッド間の間隔
w = (zmax-zmin)/(N-1);

Z = linspace(zmin,zmax,N)';
% 定常状態はmu
Z = Z + mu;

j = 4;
k = 6;

z0 = (Z(k)-c-rho*Z(j))/sigma;
z1 = (Z(k)-c-rho*Z(j)-w/2)/sigma;
z2 = (Z(k)-c-rho*Z(j)+w/2)/sigma;
f0 = normpdf(z0);
f1 = normpdf(z1);
f2 = normpdf(z2);

z = linspace(-1.0,3.0,101);
f = normpdf(z);
zz = linspace(z1,z2,101);
ff = normpdf(zz);

figure;
% pos = get(gca, 'Position');
% pos(2) = pos(2)+pos(4)/2;
% pos(4) = pos(4)/2;
% set(gca, 'Position', pos)
subplot(211);
plot(z,f,'k-','LineWidth',3);
hold on;
% h1 = area(zz,ff); % why -1???
% h1.EdgeAlpha = 0;
% h1.FaceAlpha = 0.25;
% h1.EdgeColor = [1 1 1];
% h1.FaceColor = [0 0.25 0.5];
dim = [0.55,0.75,0.10,0.10];
str = {'\Phi(x_{j}+w/2-{\rho}x_{i})', '-\Phi(x_{j}-w/2-{\rho}x_{i})'};
h4 = annotation('textbox',dim,'String',str,'LineStyle','none','FontSize',10,'FontWeight','normal');

plot([z0 z0],[0 f0],'k:');
plot([z1 z1],[0 f1],'k-');
plot([z2 z2],[0 f2],'k-');
xlabel("x'"); ylabel("\phi(x'-{\rho}x_{i})");
xticks([0 z1 z0 z2]);
xticklabels({"{\rho}x_{i}","x_{j}-w/2","x_{j}","x_{j}+w/2"});
saveas (gcf,'tauchen.eps','epsc2');