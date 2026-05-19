## 北尾早霧・砂川武貴・山田知明『定量的マクロ経済学と数値計算』日本評論社

### 第4章：オイラー方程式と時間反復法

本リポジトリは、北尾早霧・砂川武貴・山田知明『[定量的マクロ経済学と数値計算](https://www.nippyo.co.jp/shop/book/9287.html)』（日本評論社）の**第4章**に対応するサポートコードを収録しています。

---

### フォルダ構成

```
chapter4/
├── Julia/
│   ├── 4_1_4_TI_Optimization/
│   │   └── TI_Optimization.ipynb       # 4.1.4節：時間反復法（最適化）
│   ├── 4_2_1_TI_Optimization_Chebyshev/
│   │   └── TI_Optimization_Cheb.ipynb  # 4.2.1節：チェビシェフ多項式近似 + TI
│   ├── 4_2_2_TI_Optimization_EGM/
│   │   └── TI_Optimization_EGM.ipynb   # 4.2.2節：内生的格子法（EGM）+ TI
│   ├── 4_2_4_TI_StochasticOptimization/
│   │   └── TI_StochastiicOptimization.ipynb  # 4.2.4節：確率的 TI
│   └── 4_3_PerfectForesight/
│       ├── saddlepath_shooting.ipynb   # 4.3節：完全予見・サドルパス（基本）
│       ├── saddlepath_shooting_tauc.ipynb   # 4.3節：消費税ショック版
│       ├── transition2.ipynb           # 4.3節：移行ダイナミクス
│       └── transition2_tauc.ipynb      # 4.3節：移行ダイナミクス（消費税版）
├── MATLAB/
│   ├── 4_1_4_TI_Optimization/          # main_ti.m, EulerEq.m, mu_CRRA.m
│   ├── 4_2_1_TI_Optimization_Chebyshev/ # main_ti_cheb.m, polygrid.m, polybas.m 等
│   ├── 4_2_2_TI_Optimization_EGM/      # main_ti_egm.m, f.m, mu_CRRA.m 等
│   ├── 4_2_4_TI_StochasticOptimization/ # main_ti_stoch.m, EulerEq_stoch.m, tauchen.m 等
│   └── 4_3_PerfectForesight/           # saddlepath_shooting.m, transition2.m 等
└── Python/
    ├── 4_1_4_TI_Optimization/          # main_ti.py, EulerEq.py, CRRA.py + .ipynb
    ├── 4_2_1_TI_Optimization_Chebyshev/ # TI_Optimization_Chebyshev.ipynb
    ├── 4_2_2_TI_Optimization_EGM/      # TI_Optimization_EGM.ipynb
    ├── 4_2_4_TI_StochasticOptimization/ # TI_StochasticOptimization.ipynb
    └── 4_3_PerfectForesight/           # saddlepath_shooting.ipynb 等
```

---

### 各言語のコードについて

#### Julia & MATLAB
* 各節の番号に対応したフォルダにコードが格納されています。
* `4_1_4_TI_Optimization`：基本的な時間反復法の結果を再現するファイル。MATLAB 版では `main_ti.m`・`EulerEq.m`・`mu_CRRA.m` が主要ファイルです。
* `4_2_1_TI_Optimization_Chebyshev`：チェビシェフ多項式近似を TI に適用。MATLAB 版では `polygrid.m`・`polybas.m` がグリッドと基底関数の生成に使われます。
* `4_2_2_TI_Optimization_EGM`：内生的格子法（EGM）を TI に組み合わせた高速実装。
* `4_2_4_TI_StochasticOptimization`：確率的生産性ショックを含む TI。Tauchen 法による遷移確率行列の生成を含みます。
* `4_3_PerfectForesight`：完全予見モデルの移行ダイナミクスをサドルパス・シューティング法で解きます。基本ケースと消費税ショック版の両方を収録しています。

#### Python
* RAの小野泰輝さんに作成していただいたPython用再現コードです。もちろん残りうるあらゆる間違いは、すべて著者たちの責任です。

---

### Notebook の内容

#### `TI_Optimization.ipynb`　―　4.1.4節：時間反復法（最適化）

時間反復法（TI）の基本実装です。オイラー方程式を関数化し、各グリッド点で次期政策関数を所与として最適な消費（貯蓄）を数値的に求めます。VFI と比較して収束が速い点が特徴です。

- オイラー方程式の残差関数の定義と `Optim` パッケージによる最大化
- TI のループ構造と収束判定（政策関数の変化量）
- VFI との収束速度・精度の比較
- 誤差指標（オイラー方程式残差）のプロット

#### `TI_Optimization_Cheb.ipynb`　―　4.2.1節：チェビシェフ多項式近似 + TI

政策関数をチェビシェフ多項式で近似し、TI と組み合わせることで連続空間での精度向上を図ります。

- チェビシェフ格子（ゼロ点）の生成と `[-1, 1]` への変数変換
- チェビシェフ多項式基底の評価と係数の推定
- 多項式近似 TI のループ実装と収束確認
- 等間隔グリッドを用いた TI との精度・速度比較

#### `TI_Optimization_EGM.ipynb`　―　4.2.2節：内生的格子法（EGM）+ TI

EGM を TI に組み込み、最適化や求根計算なしにオイラー方程式の解を解析的に求める高速実装です。

- 次期資産グリッドを外生的に固定し、オイラー方程式の右辺から今期消費を逆算
- 予算制約を使った今期資産（内生格子）の計算
- 補間を組み合わせた政策関数の更新と収束判定
- 最適化ベース TI との速度比較

#### `TI_StochastiicOptimization.ipynb`　―　4.2.4節：確率的 TI

成長モデルに確率的な TFP ショックを加え、2次元状態空間（資本 × 生産性）で TI を実装します。

- Tauchen 法による AR(1) 過程の離散近似と遷移確率行列
- 2次元グリッド上でのオイラー方程式の期待値計算
- 確率的政策関数の収束と補間による評価
- 不確実性のない確定的モデルとの政策関数の比較

#### `saddlepath_shooting.ipynb`　―　4.3節：完全予見・サドルパス

完全予見モデルの定常状態からの逸脱と移行ダイナミクスをシューティング法で計算します。位相図に基づくサドルパスの数値的な特定が中心的なトピックです。

- 動学方程式の定式化と定常状態の計算
- シューティング法の実装：初期消費を探索して終端条件（定常状態への収束）を達成
- 位相図の描画と数値的サドルパスの重ね合わせ
- 消費税ショック版（`saddlepath_shooting_tauc.ipynb`）での政策実験

---

### 注意
* MATLABではfsolveを使っているため、インストールされているライブラリによっては動かない可能性があります。
* Juliaではいくつかのパッケージを利用しています（`Optim`、`Interpolations` 等）。
  * もし実行できない場合はREPLで`]`を押したのち、`add パッケージ名`を実行してインストールしてください。
  * 図をプロットする際に日本語が文字化けする可能性があります。数値計算そのものには無関係ですが、日本語で図を表示したい場合はフォントを追加インストールする必要があります。

---

### セットアップ

Julia と Jupyter Notebook のインストール・環境設定については、[インストールと環境構築のガイド](https://quant-macro-book.github.io/setup/) を参照してください。
