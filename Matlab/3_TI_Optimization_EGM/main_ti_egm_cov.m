clear all;

% カリブレーション
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(delta=1.0のときは解析解が存在)

% 定常状態の値
ykss = (1/beta-1+delta)/alpha;
kss = ykss^(1/(alpha-1));
yss = ykss*kss;
css = yss-delta*kss;

% EGM
% kmax = 0.5;   % 資本グリッドの最大値
% kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
kmax = 2.0*kss;   % 資本グリッドの最大値
kmin = 0.1*kss;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)

%% STEP 1(a): グリッド生成
nk   = 101; %21;    % グリッドの数
% kではなく、k'についてグリッドをとる
kpgrid = linspace(kmin, kmax, nk)';
% Carrol's (2005) change of variables
% cを求めた後、資源制約式をkについて解く必要がない
wpgrid = kpgrid.^alpha + (1-delta)*kpgrid;

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
% ここでの政策関数はc'=htilde(w')となる
cfcn0 = kpgrid;
cfcn1 = zeros(nk,1);

% 繰り返し誤差を保存する変数を設定 
dif = zeros(2,maxiter);

% 計算されたc,wを保存
cgrid = zeros(nk,1);
wgrid = zeros(nk,1);

%% STEP 4: 政策関数を繰り返し計算
while (it < maxiter && dif2 > tol)

    for i = 1:nk

        kprime = kpgrid(i);
        wprime = wpgrid(i);
        cnext = cfcn0(i); % c'=htilde(w'): wprimeにおけるcnextの値
        uc = beta*mu_CRRA(cnext,gamma)*(alpha*kprime^(alpha-1)+(1-delta));
        cons = uc^(-1/gamma); % mu_CRRAの逆関数
        % w = c+k'(非線形方程式を解く必要がない)        
        wealth = cons+kprime;        
        % 計算されたc,wを保存
        cgrid(i) = cons;
        wgrid(i) = wealth;
        
    end
    
    % 計算されたc,wから、c'=htilde(w')をアップデートする
    % extrapolation
    cfcn1 = interp1(wgrid,cgrid,wpgrid,'linear','extrap');

    % 繰り返し計算誤差を確認
    dif2 = max(abs(cfcn1-cfcn0));

    % 収束途中の繰り返し計算誤差を保存
    dif(2,it) = dif2;
    
    % 政策関数をアップデート
    cfcn0 = cfcn1;

    fprintf('iteration index: %i \n', it);
    fprintf('policy function iteration error: %1.6f \n', dif2);

    it = it + 1;

end

disp(' ');
toc;

kgrid = kpgrid; % c'=h(k') は c=h(k)と同じ

%% 最終的な政策関数が得られてから貯蓄関数を計算
pfcn0 = kgrid.^alpha + (1-delta)*kgrid - cfcn0;

%% 解析的解
p_true = beta*alpha*(kgrid.^alpha);

%%
figure;
plot(kgrid, pfcn0, '-', 'Color', 'blue', 'LineWidth', 3);
hold on;
plot(kgrid, p_true, '--', 'Color', 'red', 'LineWidth', 3);
plot(kgrid, kgrid, ':', 'Color', 'black', 'LineWidth', 2);
xlabel('今期の資本保有量：k', 'FontSize', 16);
ylabel("次期の資本保有量：k'", 'FontSize', 16);
xlim([kmin kmax]);
xticks([0.05 0.1 0.2 0.3 0.4 0.5]);
xticklabels([0.05 0.1 0.2 0.3 0.4 0.5]);
legend('近似解', '解析的解', '45度線', 'Location', 'NorthWest');
grid on;
set(gca,'FontSize', 16);
% saveas(gcf,'Fig_pti2.eps','epsc2');
