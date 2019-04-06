function m = Models(nk,kmax)

%% *** カリブレーション ***
m.beta  = 0.96; % 割引因子
m.gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
m.alpha = 0.40; % 資本分配率
m.delta = 1.00; % 固定資本減耗(0.08)
%delta = 0.08;

% *** 定常状態の値(初期値に使用) ***
m.ykss = (1/m.beta-1+m.delta)/m.alpha;
m.kss = m.ykss^(1/(m.alpha-1));
m.yss = m.ykss*m.kss;
m.css = m.yss-m.delta*m.kss;

% *** 離散化用のパラメータ ***
if (nargin>=1) 
    m.nk = nk;
else
    m.nk   = 21;    % グリッドの数
end

if (nargin>=2)
    m.kmax = kmax;
else
    m.kmax = 0.5;   % 資本グリッドの最大値
end
%kmax = 0.4;   % 資本グリッドの最大値(固定資本減耗=1.0の場合、kmax=0.5にすると59行目の初期値では解けない)
%kmax = 10.0; % 資本グリッドの最大値(固定資本減耗=0.08の場合に使用)
m.kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
%========================

%% STEP 1(a): グリッド生成

m.kgrid = linspace(m.kmin, m.kmax, m.nk)';

m.maxit = 1000;    % 繰り返し計算の最大値
m.tol  = 1.0e-5; % 許容誤差(STEP 2)