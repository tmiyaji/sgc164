# 解の追跡と分岐

## ファイル

- `pathfollowing.py`
  - 解の追跡，分岐の検出，追跡する枝の切り替えアルゴリズムの実装
  - pathfollow00
  - pathfollow
- 非線形方程式の解の追跡
  - `examples.ipynb`
    - 図5.1, 5.2, 5.3を描く．
  - `zakutsu.ipynb`
    - 梁の問題
    - 図5.4(例5.6), 図5.5(例5.7)
- 平衡点
  - `lorenz.ipynb`
    - Lorenz方程式の平衡点の追跡
    - 例5.8 図5.6
    - 非自明平衡点のHopf分岐点検出
- 写像の不動点
  - `logistic.ipynb`
    - 例5.10 ロジスティック写像の不動点の周期倍分岐
    - 図5.7
  - `leslie.ipynb`
    - 例5.11 非線型Leslieモデルの周期点の周期倍分岐
    - 図5.8
- リミットサイクル
  - `subhopf.ipynb`
    - 例5.14亜臨界Hopf分岐の標準形（柔らかい亜臨界Hopf分岐の場合）
    - Hopf分岐点を検出し，リミットサイクルの枝に切り替える．リミットサイクルはfold分岐を生じる
    - 図5.9
  - `hindmarsh.ipynb`
    - 例5.15 Hindmarsh-Roseモデル
    - 図5.10, 5.11
  - `lonrenz_cycle.ipynb`
    - 例5.16 Lorenzモデル: Hopf分岐で現れる不安定周期軌道の追跡
    - 図5.12
  - `ov.ipynb`
    - 例5.17 OVモデルにおける渋滞解の追跡
    - 図5.13

## 本文で参照していない例題

- `hindmarsh_hd.ipynb`
  - 例5.15 Hindmarsh-Roseモデル
  - リミットサイクルの追跡を有限差分法で行う
- `ChafeeInfante.ipynb`
  - Chafee-Infante 方程式の離散化
