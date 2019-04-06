function [cfcn0 dif] = nti(m)

options = optimoptions('fsolve','Display','none'); % fsolveのオプション(最適化の結果を非表示にする)

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
dif2 = 1.0;      % 政策関数の繰り返し誤差
%options.TolFun = 1.0e-10; % fsolveのオプション(最適化の許容誤差)

%% 計算開始
tic

disp(' ')
disp('-+- Solve a neoclassical growth model with time iteration -+-');
disp(' ')

%% STEP 1(b): 政策関数の初期値を当て推量

% 解析解 (for k'=g(k))
p_true = m.beta*m.alpha*(m.kgrid.^m.alpha);

% 初期値が重要
cfcn0 = m.kgrid;
%cfcn0 = m.css/m.kss*m.kgrid;
%cfcn0 = m.kgrid.^m.alpha - p_true;
%cfcn0 = css*ones(nk,1);
cfcn1 = zeros(m.nk, 1);

% 収束途中の繰り返し計算誤差を保存
% 途中経過を図示する目的なので、通常は不要(むしろ遅くなるので消すべき)
dif = zeros(2,m.maxit);

%% STEP 4: 政策関数を繰り返し計算

while it < m.maxit && dif2 > m.tol

    fprintf('iteration index: %i \n', it);
    fprintf('policy function iteration error: %e\n', dif2);

    for i = 1:m.nk

        capital = m.kgrid(i);
        wealth = capital.^m.alpha + (1.-m.delta).*capital;

        % MATLABの最適化関数(fsolve)を使って各グリッド上の政策関数の値を探す
        % 最適化の初期値は古い政策関数の値
        cons = fsolve(@EulerEq2,cfcn0(i,1),options,m,capital,cfcn0);
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