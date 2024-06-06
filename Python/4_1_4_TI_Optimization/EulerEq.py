def EulerEq(m,cons,capital,interp_c):
    """
    ----------------------------------------------
    === オイラー方程式に代入した際の誤差を返す関数 ===
    ----------------------------------------------
    <input>
    ・m: パラメータ等を格納したコンストラクタ
    ・cons: 今期の消費水準
    ・capital: 今期の資本保有量(k)
    ・interp_c: 補間した消費関数
    (※matlabコードでは関数内で消費関数の補間を行っているが、pythonコード上では
    補間した消費関数を引数で先に与える)
    <output>
    ・res: オイラー方程式に代入した際の誤差
    """
    from CRRA import CRRA, mu_CRRA
    alpha, beta, delta, kgrid = m.alpha, m.beta, m.delta, m.kgrid

    wealth = capital**alpha + (1-delta)*capital
    kprime = wealth - cons
    kprime = max(kgrid[0],kprime) #トリック:k'は正の値しか取らない

    #補間した消費関数を用いて、来期の消費水準を求める。
    cnext =  interp_c(kprime) #c'=h(k')

    #オイラー方程式誤差を求める(u'(c)=1/c を mu_CRRA関数を用いて計算している)
    #res = u'(c) - beta*u'(c')*f'(k')
    res = mu_CRRA(m,cons) - beta*mu_CRRA(m,cnext)*(alpha*kprime**(alpha-1) + 1-delta)

    return res

    

    



    