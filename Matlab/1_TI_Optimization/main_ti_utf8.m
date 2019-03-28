%% メインファイル:
% 時間反復法(time iteration method)の解法
% 状態変数のみ離散化して操作変数は連続的に値を取る(main_ndp.mと同じ)

clear;
clear global;
close all;
format short;

global beta gamma alpha delta capital cfcn kgrid
options = optimoptions('fsolve','Display','none'); % fsolveのオプション(最適化の結果を非表示にする)

%% *** カリブレーション ***
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(0.08)
%delta = 0.08;

% *** 定常状態の値(初期値に使用) ***
ykss = (1/beta-1+delta)/alpha;
kss = ykss^(1/(alpha-1));
yss = ykss*kss;
css = yss-delta*kss;

% *** 離散化用のパラメータ ***
nk   = 21;    % グリッドの数
kmax = 0.5;   % 資本グリッドの最大値
%kmax = 0.4;   % 資本グリッドの最大値(固定資本減耗=1.0の場合、kmax=0.5にすると59行目の初期値では解けない)
%kmax = 10.0; % 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
%========================

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
maxit = 1000;    % 繰り返し計算の最大値
tol  = 1.0e-5; % 許容誤差(STEP 2)
%dif1 = 1.0;      % 価値関数の繰り返し誤差
dif2 = 1.0;      % 政策関数の繰り返し誤差
%options.TolFun = 1.0e-10; % 最適化の許容誤差
count = 1;
%=================

%% 計算開始

tic

disp('')
disp('-+- Solve a neoclassical growth model with time iteration -+-');

%% STEP 1(a): グリッド生成

kgrid = linspace(kmin, kmax, nk)';

%% STEP 1(b): 政策関数の初期値を当て推量

% 初期値が重要
cfcn0 = kgrid;
p_true = beta*alpha*(kgrid.^alpha);
%cfcn0 = kgrid.^alpha - p_true;
%cfcn0 = css*ones(nk,1);
%cfcn0 = css/kss*kgrid; %zeros(nk, 1);
%vfcn0 = CRRA(kgrid.^alpha + (1.-delta).*kgrid, gamma);

cfcn1 = zeros(nk, 1);
%vfcn1 = zeros(nk, 1);

% 価値関数・政策関数の経路を記録(なくても可)
%vpath(:, 1) = vfcn0;
ppath(:, 1) = cfcn0;

%% STEP 4: 政策関数を繰り返し計算

while it < maxit && dif2 > tol

    fprintf('iteration index: %i \n', it);
%    fprintf('value function iteration error: %e\n', dif1);
    fprintf('policy function iteration error: %e\n', dif2);

    for i = 1:nk

        % グローバル変数を設定
        % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
        capital = kgrid(i);
        wealth = capital.^alpha + (1.-delta).*capital;
        cfcn = cfcn0;

        % MATLABの最適化関数(fsolve)を使ってグリッド上で政策関数の値を探す
        % 最適化の初期値は古い政策関数の値
        cons = fsolve(@EulerEq, cfcn0(i,1), options);
        % 最適化の初期値は定常状態の値: これでは解けない
        % cons = fsolve(@EulerEq, css, options);
        cfcn1(i,1) = cons;
        kprime = wealth-cons;
%        disp([cons capital wealth kprime]);
%        pause

    end

    % fminsearchは最小値を探す関数なので符号を反転させる
%    vfcn1 = -1*vfcn1;

    % 繰り返し計算誤差を確認
%    dif1 = max(abs((vfcn1-vfcn0)./vfcn0));
%    dif2 = max(abs((cfcn1-cfcn0)./cfcn0));
    dif2 = max(abs(cfcn1-cfcn0));

    % 収束途中の繰り返し計算誤差を保存
    % 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
    % 計算毎に行列のサイズが変わっていくのは望ましくない書き方なので本来は避けるべき
%    dif(1, it) = dif1;
    dif(2, it) = dif2;

    % 価値関数の経路を記録(なくても問題なし)
%    vpath(:, it) = vfcn0;
    ppath(:, it) = cfcn0;

    % 価値関数・政策関数をアップデート
%    vfcn0 = vfcn1;
    cfcn0 = cfcn1;

    it = it + 1;

end

% 最終的な政策関数が得られてから貯蓄関数を計算
pfcn0 = kgrid.^alpha + (1-delta)*kgrid - cfcn0;

%% 計算結果をコマンドウィンドウに表示

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

toc

%% 解析的解

AA = (1.0-beta).^(-1) * (log(1.0-alpha*beta) + ((alpha*beta)/(1.0-alpha*beta))*log(alpha*beta));
BB = alpha/(1.0-alpha*beta);
v_true = AA + BB*log(kgrid);
p_true = beta*alpha*(kgrid.^alpha);

%% オイラー方程式から誤差を測定
% 元のグリッドではオイラー方程式の誤差はゼロになるため、グリッドを細かくとる
kgrid_err = linspace(kmin, kmax, (nk-1)*10+1)';
%kgrid_err = kgrid;
%cons = kgrid.^alpha + (1.-delta).*kgrid - pfcn0(:,1);
cons = interp1(kgrid,cfcn0(:,1),kgrid_err);
LHS  = mu_CRRA(cons, gamma);
kp   = kgrid_err.^alpha + (1-delta)*kgrid_err - cons;
cnext = interp1(kgrid, cfcn0(:,1), kp);
%cnext = kp.^alpha + (1.-delta).*kp - kpp;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cnext, gamma);
err  = RHS./LHS-1.0;
% 
% err2 = csvread("err_ddp.csv");
% 
%% 図を描く
% 
% figure;
% plot(kgrid, vfcn0, '-', 'color', 'blue', 'linewidth', 3); hold('on');
% plot(kgrid, v_true, '--', 'color', 'red', 'linewidth', 3); hold('off');
% %title('価値関数', 'fontsize', 16);
% xlabel('資本保有量：k', 'Fontsize', 16);
% ylabel('価値関数：V(k)', 'Fontsize', 16);
% xlim([0, kmax]);
% legend('近似解', '解析的解', 'Location', 'SouthEast');
% grid on;
% set(gca,'Fontsize', 16);
% saveas (gcf,'Fig3_pndp1.eps','epsc2');
% 
figure;
plot(kgrid, pfcn0, '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kgrid, p_true, '--', 'color', 'red', 'linewidth', 3);
plot(kgrid, kgrid, ':', 'color', 'black', 'linewidth', 2); hold('off');
%title('政策関数', 'fontsize', 16);
xlabel('今期の資本保有量：k', 'Fontsize', 16);
ylabel("次期の資本保有量：k'", 'Fontsize', 16);
xlim([0, kmax]);
legend('近似解', '解析的解', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig_pti2.eps','epsc2');
 
figure;
plot(kgrid, cfcn(:,1), '-', 'color', 'blue', 'linewidth', 3);
%title('消費関数', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('消費：c', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig_pti3.eps','epsc2');
 
iter = linspace(1, it-1, it-1);
 
figure;
% plot(iter, dif(1, :), '-', 'color', 'blue', 'linewidth', 2); hold('on');
plot(iter, dif(2, :), ':', 'color', 'red', 'linewidth', 2); hold('off');
% %title('価値関数・政策関数の収束', 'fontsize', 16);
xlabel('計算回数', 'Fontsize', 16);
ylabel('繰り返し計算誤差', 'Fontsize', 16);
ylim([0,0.1]);
% legend('価値関数', '政策関数', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig_pti4.eps','epsc2');

figure;
plot(kgrid_err, err, '-', 'color', 'blue', 'linewidth', 3);
%title('オイラー方程式の誤差', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('オイラー方程式誤差', 'Fontsize', 16);
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig_pti5.eps','epsc2');

err2 = csvread("err_ndp.csv");
figure;
plot(kgrid_err, abs(err), '-', 'color', 'blue', 'linewidth', 3); hold('on');
plot(kgrid_err, abs(err2), '--', 'color', 'red', 'linewidth', 3); hold('off');
%title('オイラー方程式の誤差', 'fontsize', 16);
xlabel('資本保有量：k', 'Fontsize', 16);
ylabel('オイラー方程式誤差(絶対値)', 'Fontsize', 16);
%ylim([-15e-004,5e-004]);
legend('TI', 'VFI', 'Location', 'NorthEast');
grid on;
set(gca,'Fontsize', 16);
saveas (gcf,'Fig_pti6.eps','epsc2');

return
