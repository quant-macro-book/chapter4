clear all;

% MATLABではグローバル変数を使って関数への変数受け渡しを行う
% Julia, Pythondではコンストラクタを用いる
global beta gamma alpha delta kgrid zgrid Pzvec

% カリブレーション
beta  = 0.99; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.36; %0.40; % 資本分配率
delta = 0.025; %1.00; % 固定資本減耗(delta=1.0のときは解析解が存在)

% 定常状態の値
ykss = (1/beta-1+delta)/alpha;
kss = ykss^(1/(alpha-1));
yss = ykss*kss;
css = yss-delta*kss;

%kmax = 0.5;   % 資本グリッドの最大値
%kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
% kmax = 1.2*kss;   % 資本グリッドの最大値
% kmin = 0.8*kss;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
kmax = 1.5*kss;   % 資本グリッドの最大値
kmin = 0.1*kss;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)

%% STEP 1(a): グリッド生成
nk   = 21;    % グリッドの数
kgrid = linspace(kmin, kmax, nk)';
nz   = 11;
rho = 0.95;
sigma = 0.01;
m = 2.575;
[zgrid,Pz] = tauchen(nz,0,rho,sigma,m);

maxiter = 1000; % 繰り返し計算の最大値
tol  = 1.0e-5;  % 許容誤差(STEP 2)

% 収束の基準に関するパラメータ
it = 1;          % ループ・カウンター
dif2 = 1.0;      % 政策関数の繰り返し誤差
options.TolFun = 1.0e-10; % fsolveのオプション(最適化の許容誤差)

tic;
options = optimoptions('fsolve','Display','none'); % fsolveのオプション(最適化の結果を非表示にする)

disp(' ')
disp('-+- Solve a neoclassical growth model with time iteration -+-');
disp(' ')

%% STEP 1(b): 政策関数の初期値を当て推量
cfcn0 = (css/kss)*kgrid*ones(1,nz);
cfcn1 = zeros(nk,nz);
mfcn0 = mu_CRRA(css,gamma)*(alpha*kgrid.^(alpha-1) + (1.-delta))*ones(1,nz);
mfcn1 = zeros(nk,nz);

% 繰り返し誤差を保存する変数を設定 
dif = zeros(2,maxiter);

%% STEP 4: 政策関数を繰り返し計算
while (it < maxiter && dif2 > tol)

    for iz = 1:nz
        
        technology = exp(zgrid(iz));
        Pzvec = Pz(iz,:)';
        mfcn0 = mfcn0(:,iz)*Pzvec;
        
        for ik = 1:nk

            capital = kgrid(ik);
            wealth = technology*capital.^alpha + (1.-delta).*capital;

            % MATLABの最適化関数(fsolve)を使って各グリッド上の政策関数の値を探す
            % 最適化の初期値は古い政策関数の値
            cons = fsolve(@EulerEq_stoch,cfcn0(ik,iz),options,wealth,cfcn0);
            cfcn1(ik,iz) = cons;
            kprime = wealth-cons;
            % change of variables
            mfcn1(ik,iz) = mu_CRRA(cons,gamma)*(alpha*exp(zgrid(iz))*capital.^(alpha-1) + (1.-delta));
            % グリッドごとに最適化の結果を確認
            %disp([cons capital wealth kprime]);
            %pause

        end
    
    end

    % 繰り返し計算誤差を確認
    dif2 = max(max(abs(cfcn1-cfcn0)));

    % 収束途中の繰り返し計算誤差を保存
    dif(2,it) = dif2;
    
    % 政策関数をアップデート
    cfcn0 = cfcn1;
%    cfcn0 = 0.1*cfcn1 + 0.9*cfcn0;

    fprintf('iteration index: %i \n', it);
    fprintf('policy function iteration error: %1.6f \n', dif2);

    it = it + 1;

end

disp(' ');
toc;

% consumption path
% simulation
% generate exogenous shock sequence
T = 1000;
ivec = zeros(T+1,1);
ivec(1) = 3;
% 条件付き累積密度関数：i列の要素の累積和
cumP = cumsum(Pz')';

for t = 1:T
    cumPi = cumP(ivec(t),:);
    % 一様分布から[0,1]の値をとる乱数を発生させ、条件付き累積密度関数と比較する
    % たとえば、rand<cumPi(1)のとき、ivec(t+1)=1となる
    ivec(t+1) = sum(rand-cumPi >= 0)+1;
    ivec(t+1) = min(ivec(t+1),nz);
end

cvec = zeros(T,1);
kvec = zeros(T+1,1);
kvec(1) = kss;
for t = 1:T
    iz = ivec(t);
    know = kvec(t);
    cnow = interp1(kgrid,cfcn0(:,iz),know,'linear','extrap');
    kprime = exp(zgrid(iz))*know^alpha + (1-delta)*know - cnow;
    cvec(t) = cnow;
    kvec(t+1) = kprime;
end

figure;
subplot(311);
plot(cvec);
subplot(312);
plot(kvec(1:T));
subplot(313);
plot(exp(zgrid(ivec(1:T))));

% %% 最終的な政策関数が得られてから貯蓄関数を計算
% pfcn0 = kgrid.^alpha + (1-delta)*kgrid - cfcn0;
% 
% %% 解析的解
% p_true = beta*alpha*(kgrid.^alpha);
% 
% %%
% figure;
% plot(kgrid, pfcn0, '-', 'Color', 'blue', 'LineWidth', 3);
% hold on;
% plot(kgrid, p_true, '--', 'Color', 'red', 'LineWidth', 3);
% plot(kgrid, kgrid, ':', 'Color', 'black', 'LineWidth', 2);
% xlabel('今期の資本保有量：k', 'FontSize', 16);
% ylabel("次期の資本保有量：k'", 'FontSize', 16);
% xlim([kmin kmax]);
% xticks([0.05 0.1 0.2 0.3 0.4 0.5]);
% xticklabels([0.05 0.1 0.2 0.3 0.4 0.5]);
% legend('近似解', '解析的解', '45度線', 'Location', 'NorthWest');
% grid on;
% set(gca,'FontSize', 16);
% saveas(gcf,'Fig_pti2.eps','epsc2');
% 
% %% オイラー方程式から誤差を測定
% % 元のグリッドではオイラー方程式の誤差はゼロになるため、グリッドを細かくとる
% kgrid_err = linspace(kmin, kmax, (nk-1)*10+1)';
% % c=h(k;b)の値を線形補間
% cons = interp1(kgrid,cfcn0(:,1),kgrid_err);
% LHS  = mu_CRRA(cons, gamma);
% 
% % k'=f(k)-c
% kp   = kgrid_err.^alpha + (1-delta)*kgrid_err - cons;
% % c'=h(k';b)の値を線形補間
% cnext = interp1(kgrid, cfcn0(:,1), kp);
% rent = alpha.*kp.^(alpha-1.0) - delta;
% RHS  = beta.*(1.+rent).*mu_CRRA(cnext,gamma);
% 
% err  = RHS./LHS-1.0;
% 
% err2 = csvread("err_ndp.csv");
% figure;
% plot(kgrid_err, abs(err), '-', 'Color', 'blue', 'LineWidth', 3);
% hold on;
% plot(kgrid_err, abs(err2), '--', 'Color', 'red', 'LineWidth', 3);
% xlabel('資本保有量：k', 'FontSize', 16);
% ylabel('オイラー方程式誤差(絶対値)', 'FontSize', 16);
% xlim([kmin kmax]);
% xticks([0.05 0.1 0.2 0.3 0.4 0.5]);
% xticklabels([0.05 0.1 0.2 0.3 0.4 0.5]);
% legend('TI', 'VFI', 'Location', 'NorthEast');
% grid on;
% set(gca,'FontSize', 16);
% saveas (gcf,'Fig_pti6.eps','epsc2');
