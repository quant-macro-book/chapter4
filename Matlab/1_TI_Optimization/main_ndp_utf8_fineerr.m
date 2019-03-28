%% メインファイル:
% 状態変数のみ離散化して操作変数は連続的に値を取る場合の動的計画法(parametric DP)の解法.
% アルゴリズムの詳細は、Johnson et al. (1993)を参照

clear;
clear global;
close all;
format short;

global beta gamma alpha delta capital vfcn kgrid

%% *** カリブレーション ***
beta  = 0.96; % 割引因子
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(0.08)

% *** 離散化用のパラメータ ***
nk   = 21;    % グリッドの数
kmax = 0.5;   % 資本グリッドの最大値
%kmax = 10.0; % 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
%========================

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
maxit = 1000;    % 繰り返し計算の最大値
tol  = 1.0e-005; % 許容誤差(STEP 2)
dif1 = 1.0;      % 価値関数の繰り返し誤差
dif2 = 1.0;      % 政策関数の繰り返し誤差
count = 1;
%=================

%% 計算開始

tic

disp('')
disp('-+- Solve a neoclassical growth model -+-');

%% STEP 1(a): グリッド生成

kgrid = linspace(kmin, kmax, nk)';
% kgrid = grid_exp1(kmin, kmax, nk)';
% kgrid = grid_exp2(kmin, kmax, nk)';
% kgrid = grid_exp3(kmin, kmax, nk)';

%% STEP 1(b): 価値関数・政策関数の初期値を当て推量

pfcn0 = zeros(nk, 1);
vfcn0 = CRRA(kgrid.^alpha + (1.-delta).*kgrid, gamma);

pfcn1 = zeros(nk, 1);
vfcn1 = zeros(nk, 1);

% 価値関数・政策関数の経路を記録(なくても可)
vpath(:, 1) = vfcn0;
ppath(:, 1) = pfcn0;

%% STEP 4: 価値関数を繰り返し計算

while it < maxit && dif1 > tol

    fprintf('iteration index: %i \n', it);
    fprintf('value function iteration error: %e\n', dif1);
    fprintf('policy function iteration error: %e\n', dif2);

    for i = 1:nk

        % グローバル変数を設定
        % fminsearchで使う関数(BellmanEq)に最適化する変数"以外"の変数を渡す(グローバル変数を使わない方法もあるはず)
        capital = kgrid(i);
        vfcn = vfcn0;

        % MATLABの最適化関数(fminsearch)を使ってグリッド上で価値関数と政策関数の値を探す
        % 初期値は0.01
        [pfcn1(i,1), vfcn1(i,1)] = fminsearch(@BellmanEq, 0.01);

    end

    % fminsearchは最小値を探す関数なので符号を反転させる
    vfcn1 = -1*vfcn1;

    % 繰り返し計算誤差を確認
    dif1 = max(abs((vfcn1-vfcn0)./vfcn0));
    dif2 = max(abs((pfcn1-pfcn0)./pfcn0));

    % 収束途中の繰り返し計算誤差を保存
    % 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
    % 計算毎に行列のサイズが変わっていくのは望ましくない書き方なので本来は避けるべき
    dif(1, it) = dif1;
    dif(2, it) = dif2;

    % 価値関数の経路を記録(なくても問題なし)
    vpath(:, it) = vfcn0;
    ppath(:, it) = pfcn0;

    % 価値関数・政策関数をアップデート
    vfcn0 = vfcn1;
    pfcn0 = pfcn1;

    it = it + 1;

end

% 最終的な政策関数が得られてから消費関数を計算
cfcn = kgrid.^alpha + (1.-delta).*kgrid - pfcn0(:,1);

%% 計算結果をコマンドウィンドウに表示

disp('-+- PARAMETER VALUES -+-');
disp('');
fprintf('beta=%5.2f, gamma=%5.2f, alpha=%5.2f, delta=%5.2f \n', beta, gamma, alpha, delta);
disp(''); 
fprintf('kmin=%5.2f, kmax=%5.2f, #grid=%i \n', kmin, kmax, nk);
disp('');

toc

%% オイラー方程式から誤差を測定
% 元のグリッドではオイラー方程式の誤差はゼロになるため、グリッドを細かくとる
kgrid_err = linspace(kmin, kmax, (nk-1)*10+1)';
kp   = interp1(kgrid,pfcn0(:,1),kgrid_err);
cons = kgrid_err.^alpha + (1.-delta).*kgrid_err - kp;
LHS  = mu_CRRA(cons, gamma);
kpp  = interp1(kgrid, pfcn0(:,1), kp);
cons = kp.^alpha + (1.-delta).*kp - kpp;
rent = alpha.*kp.^(alpha-1.0) - delta;
RHS  = beta.*(1.+rent).*mu_CRRA(cons, gamma);
err  = RHS./LHS-1.0;

csvwrite("err_ndp.csv",err);

return
