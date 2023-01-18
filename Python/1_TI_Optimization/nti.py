def nti(m):
    """
    ---------------------------
    === 時間反復法を行う関数 ===
    ---------------------------
    <input>
    ・m: パラメータ等を格納したコンストラクタ
    <output>
    ・cfcn0: 収束した消費関数
    """
    import numpy as np
    from scipy.optimize import fsolve, newton
    from scipy.interpolate import interp1d
    from EulerEq import EulerEq

    kgrid, nk, alpha, delta = m.kgrid, m.nk, m.alpha, m.delta
    tol, maxiter = m.tol, m.maxiter

    #収束の基準に関するパラメータ
    it = 1    #ループ・カウンター
    dif = 1.0 #政策関数の繰り返し誤差

    #STEP 1(b): 政策関数の初期値を当て推量 
    cfcn0 = kgrid
    cfcn1 = np.zeros(nk)

    #STEP 4: 政策関数を繰り返し計算
    while (it < maxiter) & (dif > tol):

        interp_c = interp1d(kgrid,cfcn0,kind="cubic",fill_value="extrapolate")

        for i in range(nk):

            capital = kgrid[i]
            wealth = (capital ** alpha) + (1-delta)*capital

            #Pythonの最適化関数(newton)を使って各グリッド上の政策関数の値を探す
            #-> オイラー方程式の誤差をゼロにするようなcの値を求める
            #最適化の初期値は古い政策関数の値
            Euler = lambda x: EulerEq(m,x,capital,interp_c)
            cfcn1[i] = newton(Euler, x0 = cfcn0[i], tol=1e-10)

        # 繰り返し計算誤差を確認
        dif = np.max(np.abs(cfcn1-cfcn0))

        #政策関数をアップデート
        cfcn0 = np.copy(cfcn1)

        #繰り返しの結果を表示
        print(f"iteration index: {it}")
        print(f"policy function iteration error: {dif}")

        it += 1 


    return cfcn0     



