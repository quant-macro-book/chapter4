clear all;

% カリブレーション
m.beta  = 0.96; % 割引因子
m.gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
m.alpha = 0.40; % 資本分配率
m.delta = 1.00; % 固定資本減耗(delta=1.0のときは解析解が存在)

% 定常状態の値
m.ykss = (1/m.beta-1+m.delta)/m.alpha;
m.kss = m.ykss^(1/(m.alpha-1));
m.yss = m.ykss*m.kss;
m.css = m.yss-m.delta*m.kss;

m.kmax = 0.5;   % 資本グリッドの最大値
m.kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)

%% STEP 1(a): グリッド生成
m.nk   = 21;    % グリッドの数
m.kgrid = linspace(m.kmin, m.kmax, m.nk)';

m.maxiter = 1000; % 繰り返し計算の最大値
m.tol  = 1.0e-5;  % 許容誤差(STEP 2)

tic;

options = optimoptions('fsolve','Display','none'); % fsolveのオプション(最適化の結果を非表示にする)

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
dif2 = 1.0;      % 政策関数の繰り返し誤差
%options.TolFun = 1.0e-10; % fsolveのオプション(最適化の許容誤差)

disp(' ')
disp('-+- Solve a neoclassical growth model with time iteration -+-');
disp(' ')

%% STEP 1(b): 政策関数の初期値を当て推量
% 解析解 (for k'=g(k))
p_true = m.beta*m.alpha*(m.kgrid.^m.alpha);

% 政策関数の初期化
cfcn0 = m.kgrid;
%cfcn0 = m.css/m.kss*m.kgrid; % m.nk=21のときは政策関数の形がおかしい???
%cfcn0 = m.kgrid.^m.alpha - p_true;
%cfcn0 = m.css*ones(nk,1);
cfcn1 = zeros(m.nk,1);

% 繰り返し誤差を保存する変数を設定 
dif = zeros(2,m.maxiter);

%% STEP 4: 政策関数を繰り返し計算
while (it < m.maxiter && dif2 > m.tol)

    fprintf('iteration index: %i \n', it);
    fprintf('policy function iteration error: %e\n', dif2);

    for i = 1:m.nk

        capital = m.kgrid(i);
        wealth = capital.^m.alpha + (1.-m.delta).*capital;

        % MATLABの最適化関数(fsolve)を使って各グリッド上の政策関数の値を探す
        % 最適化の初期値は古い政策関数の値
        cons = fsolve(@EulerEq,cfcn0(i,1),options,m,capital,cfcn0);
        % 最適化の初期値は定常状態の値: これでは解けない
        % cons = fsolve(@EulerEq2,css,options,m,capital,cfcn0);
        cfcn1(i,1) = cons;
        kprime = wealth-cons;
        % グリッドごとに最適化の結果を確認
        %disp([cons capital wealth kprime]);
        %pause

    end

    % 繰り返し計算誤差を確認
    dif2 = max(abs(cfcn1-cfcn0));

    % 収束途中の繰り返し計算誤差を保存
    dif(2,it) = dif2;

    % 政策関数をアップデート
    cfcn0 = cfcn1;

    it = it + 1;

end

disp(' ');
toc;

%% 最終的な政策関数が得られてから貯蓄関数を計算
pfcn0 = m.kgrid.^m.alpha + (1-m.delta)*m.kgrid - cfcn0;

%% 解析的解
p_true = m.beta*m.alpha*(m.kgrid.^m.alpha);

%% オイラー方程式から誤差を測定
% 元のグリッドではオイラー方程式の誤差はゼロになるため、グリッドを細かくとる
kgrid_err = linspace(m.kmin, m.kmax, (m.nk-1)*10+1)';
cons = interp1(m.kgrid,cfcn0(:,1),kgrid_err); % 線形補間
LHS  = mu_CRRA(cons, m.gamma);

kp   = kgrid_err.^m.alpha + (1-m.delta)*kgrid_err - cons;
cnext = interp1(m.kgrid, cfcn0(:,1), kp);
rent = m.alpha.*kp.^(m.alpha-1.0) - m.delta;
RHS  = m.beta.*(1.+rent).*mu_CRRA(cnext,m.gamma);

err  = RHS./LHS-1.0;

%%
figure;
plot(m.kgrid, pfcn0, '-', 'Color', 'blue', 'LineWidth', 3);
hold on;
plot(m.kgrid, p_true, '--', 'Color', 'red', 'LineWidth', 3);
plot(m.kgrid, m.kgrid, ':', 'Color', 'black', 'LineWidth', 2);
xlabel('今期の資本保有量：k', 'FontSize', 16);
ylabel("次期の資本保有量：k'", 'FontSize', 16);
xlim([m.kmin m.kmax]);
xticks([0.05 0.1 0.2 0.3 0.4 0.5]);
xticklabels([0.05 0.1 0.2 0.3 0.4 0.5]);
legend('近似解', '解析的解', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'FontSize', 16);
saveas(gcf,'Fig_pti2.eps','epsc2');

err2 = csvread("err_ndp.csv");
figure;
plot(kgrid_err, abs(err), '-', 'Color', 'blue', 'LineWidth', 3);
hold on;
plot(kgrid_err, abs(err2), '--', 'Color', 'red', 'LineWidth', 3);
xlabel('資本保有量：k', 'FontSize', 16);
ylabel('オイラー方程式誤差(絶対値)', 'FontSize', 16);
xlim([m.kmin m.kmax]);
xticks([0.05 0.1 0.2 0.3 0.4 0.5]);
xticklabels([0.05 0.1 0.2 0.3 0.4 0.5]);
legend('TI', 'VFI', 'Location', 'NorthEast');
grid on;
set(gca,'FontSize', 16);
saveas (gcf,'Fig_pti6.eps','epsc2');
