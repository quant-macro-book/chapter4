#定量的マクロ経済学と数値計算 第4章 オイラー方程式と多項式近似
#時間反復法: メインファイル

#必要な関数・モジュールを呼び出す
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import japanize_matplotlib
from CRRA import CRRA, mu_CRRA
from EulerEq import EulerEq
from nti import nti

#カリブレーション
class Model():
    """
    モデルを解くためのパラメータを含む class を定義する。
    """
    def __init__(self,
        beta = 0.96,    # 割引因子
        gamma = 1.0,    # 相対的リスク回避度(異時点間の代替弾力性の逆数)
        alpha = 0.4,    # 資本分配率 
        delta = 1.0,    # 固定資本減耗
        nk = 21,        # 資本のグリッド数
        kmax = 0.5,     # 資本グリッドの最大値
        kmin = 0.05,    # 資本グリッドの最小値
        maxiter = 1000, # 繰り返し計算の最大値
        tol = 1e-5,     # 許容誤差(STEP2)
        ): 
        
        self.beta, self.gamma, self.alpha = beta, gamma, alpha 
        self.delta, self.nk = delta, nk 
        self.kmax, self.kmin = kmax, kmin 
        self.kgrid = np.linspace(kmin,kmax,nk) # 資本のグリッド
        self.maxiter, self.tol = maxiter, tol

print("")
print("-+- Solve a neoclassical growth model with time iteration -+-")
print("")

m = Model()
mu_CRRA_vec = np.vectorize(mu_CRRA)

#時間反復法を実施
cfcn0 = nti(m)

#最終的な政策関数が得られてから貯蓄関数を計算
pfcn0 =  (m.kgrid ** m.alpha) + (1-m.delta)*m.kgrid - cfcn0 

#解析的解(for k'=g(k))
p_true = m.beta * m.alpha * (m.kgrid ** m.alpha)

#オイラー方程式から誤差を測定
kgrid_err = np.linspace(m.kmin,m.kmax,(m.nk-1)*10+1)
cons_interp = interp1d(m.kgrid,cfcn0,kind="linear",fill_value="extrapolate")
cons = cons_interp(kgrid_err)
LHS = mu_CRRA_vec(m,cons)

kp = (kgrid_err ** m.alpha) + (1-m.delta)*kgrid_err - cons
cnext = cons_interp(kp)
rent = m.alpha*(kp ** (m.alpha-1)) - m.delta
RHS = m.beta * (1+rent) * mu_CRRA_vec(m,cnext)

err = (RHS/LHS) - 1

#図を描く
fig, ax = plt.subplots()
ax.plot(m.kgrid,pfcn0,label="近似解")
ax.plot(m.kgrid,p_true,ls="--",label="解析的解")
ax.plot(m.kgrid,m.kgrid,label="45度線")
ax.set(xlabel=r"今期の資本保有量: $k$",ylabel=r"次期の資本保有量: $k'$",title=r"(a):政策関数:$g(x)=f(x)-h(x)$",
xlim=(m.kmin,m.kmax),xticks=[0.05,0.1,0.2,0.3,0.4,0.5])
ax.legend(loc="upper left")
ax.grid()
plt.show()

f = open("err_ndp.csv",encoding='utf-8-sig')
err2 = np.loadtxt(f) #VFIでのオイラー方程式誤差

fig, ax = plt.subplots()
ax.plot(kgrid_err,np.abs(err),label="TI")
ax.plot(kgrid_err,np.abs(err2),ls="--",label="VFI")
ax.set(xlabel=r"資本保有量: $k$",ylabel="オイラー方程式誤差(絶対値)",title="(b):オイラー方程式に代入した際の誤差",
xlim=(m.kmin,m.kmax),xticks=[0.05,0.1,0.2,0.3,0.4,0.5])
ax.legend(loc="upper right")
ax.grid()
plt.show()

#TIとVFIでのオイラー方程式誤差を比較する
L1_TI = np.log10(np.mean(np.abs(err)))
Lmax_TI = np.log10(np.max(np.abs(err)))

L1_VFI = np.log10(np.mean(np.abs(err2)))
Lmax_VFI = np.log10(np.max(np.abs(err2)))

print("")
print("method", "[log10(L1), log10(L∞)]")
print("TI", [L1_TI, Lmax_TI])
print("VFI", [L1_VFI, Lmax_VFI])
print("")
