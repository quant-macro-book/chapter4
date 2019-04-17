%clear all;

%%
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

% m.kmax = 0.5;   % 資本グリッドの最大値
% m.kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
m.kmax = 1.2*m.kss;  % 資本グリッドの最大値
m.kmin = 0.8*m.kss;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)

m.maxiter = 1000; % 繰り返し計算の最大値
m.tol  = 1.0e-8;  % 許容誤差(STEP 2)

%%
norms = zeros(3,2);
times = zeros(3,2);

nkvec = [3 5 9]';

for i=1:3

    %% STEP 1(a): グリッド生成
    m.nk = nkvec(i);
    m.kgrid = polygrid(m.kmin,m.kmax,m.nk);
    m.T = polybas(m.kmin,m.kmax,m.nk,m.kgrid);
    m.invT = inv(m.T);

    % time iteration
    tic;
    [cfcn0 dif] = nti_cheb(m);
    times(i,1) = toc;

    err = calcerr(m,cfcn0);
    norms(i,:) = log10([mean(abs(err)) max(abs(err))]);
    
end

times(:,2) = times(:,1)/times(1,1);

disp(" Euler equation errors");
disp([round(norms,2)]);
disp(" Elasped time");
disp([round(times,2)]);