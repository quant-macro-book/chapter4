clear all;

% MATLABではグローバル変数を使って関数への変数受け渡しを行う
% Julia, Pythonではコンストラクタを用いる
global beta gamma alpha delta kmin kmax nk kgrid

% カリブレーション
beta  = 0.96; % 割引因子
%beta = 0.99;
gamma = 1.0;  % 相対的危険回避度(異時点間の代替の弾力性の逆数)
alpha = 0.40; % 資本分配率
delta = 1.00; % 固定資本減耗(delta=1.0のときは解析解が存在)
%delta = 0.025;

% 定常状態の値
ykss = (1/beta-1+delta)/alpha;
kss = ykss^(1/(alpha-1));
yss = ykss*kss;
css = yss-delta*kss;

% kmax=0.5のときは収束しない
% kmax = 0.5;   % 資本グリッドの最大値
% kmin = 0.05;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)
kmax = 1.2*kss;  % 資本グリッドの最大値
kmin = 0.8*kss;  % 資本グリッドの最小値 (0にすると生産が出来なくなる)

%% STEP 1(a): グリッド生成
nk = 5;    % グリッドの数
% チェビシェフ極値点あるいはゼロ点をグリッドとする
kgrid = polygrid(kmin,kmax,nk);
% グリッド上で基底関数を評価し、その逆行列をとる
T = polybas(kmin,kmax,nk,kgrid);
invT = inv(T);

maxiter = 1000; % 繰り返し計算の最大値
tol  = 1.0e-8;  % 許容誤差(STEP 2)

% *** 収束の基準 ***
it = 1;          % ループ・カウンター
dif2 = 1.0;      % 政策関数の繰り返し誤差
options.TolFun = 1.0e-10; % fsolveのオプション(最適化の許容誤差)

tic;
options = optimoptions('fsolve','Display','none'); % fsolveのオプション(最適化の結果を非表示にする)

disp(' ')
disp('-+- Solve a neoclassical growth model with projection method -+-');
disp(' ')

%% STEP 1(b): 政策関数の初期値を当て推量
% 政策関数の初期化
cfcn0 = css/kss*kgrid;
cfcn1 = zeros(nk,1);

% fsolveの設定：(i)レーベンバーグ・マルカート法を使う、(ii)反復の最大回数を1000回に設定
options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 10000); %,'display','off');
%options = optimoptions('fminunc', 'MaxFunctionEvaluations', 100000, 'MaxIterations', 10000); %,'display','off');
damp = 1.0;

%% STEP 4: 政策関数を計算
% 補間はあらかじめグリッド上で計算した基底関数invTにより行う
% thetaは多項式の係数
theta = invT*cfcn0;

% fsolveを使って、選点上で残差がゼロに近くなる係数thetaを探す
% 注意：Student editionではfsolveが入っていないかも知れません.
theta_new = fsolve(@EulerEq_cheb_proj_noloop,theta,options,T);
%theta_new = fminunc(@EulerEq_cheb_proj_noloop2,theta,options,T);
cfcn1 = T*theta_new;

% 政策関数をアップデート
cfcn0 = damp*cfcn1 + (1-damp)*cfcn0;
theta = damp*theta_new + (1-damp)*theta;

disp(' ');
toc;

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

err = calcerr(invT,cfcn0);
disp(" Euler equation errors");
disp(log10([mean(abs(err)) max(abs(err))]))