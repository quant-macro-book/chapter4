def CRRA(m,cons):
        """
        ------------------------------
        === CRRA Utility Function ===
        ------------------------------
        <inputs>
        ・m: パラメータ等を格納したコンストラクタ
        ・cons: 消費量
        ・gamma: 異時点間の代替弾力性の逆数(相対的リスク回避度)
        <output>
        ・consとgamma の下での効用水準
        """
        import numpy as np
        gamma = m.gamma
        
        if cons <0: 
            utility = -10000 #消費量が負値の場合、ペナルティを与える(最適化問題におけるトリック)
        else:
            if gamma != 1:
                utility = (cons ** (1-gamma)) / (1-gamma)
            else:
                utility = np.log(cons)
                
        return utility

def mu_CRRA(m,cons):
        """
        --------------------------------------
        === CRRA Marginal Utility Function ===
        --------------------------------------
        <inputs>
        ・m: パラメータ等を格納したコンストラクタ
        ・cons: 消費量
        <output>
        ・consとgamma の下での限界効用水準
        """
        gamma = m.gamma
        mu = cons ** (-gamma)

        return mu
