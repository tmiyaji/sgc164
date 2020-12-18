# 「数理モデルとシミュレーション」サポートページ

SGCライブラリ  164 「数理モデルとシミュレーション」
小川知之，宮路智行（共著）
サイエンス社，ISBN：978-4-7819-1495-4

## contents

- files: 書籍の本文掲載の.pyファイル
- notebooks: 本文掲載とは限らない，より詳細なコードのJupyter notebook
  - simulation
    - 「第2章 微分方程式のシミュレーション」
  - genmodel
    - 「第3章 現象論的モデル」
  - oscillation
    - 「第5章 リミットサイクル」
  - continuation
    - 「第6章 解の追跡と分岐」
  - random
    - 「第7章 擬似乱数生成器とモンテカルロ積分」

## Environment

著者が動作を確認した環境（動作を保証するものではありません）
- Ubuntu 20.04 on WSL2 on Windows 10 version 2004
  - Python 3.8.2
  - numpy 1.19.2
  - scipy 1.5.2
  - matplotlib 3.3.1
  - seaborn 0.11.0
  - `jupyter --version`
  ~~~
    jupyter core     : 4.6.3
    jupyter-notebook : 6.1.4
    qtconsole        : 4.7.7
    ipython          : 7.18.1
    ipykernel        : 5.3.4
    jupyter client   : 6.1.7
    jupyter lab      : not installed
    nbconvert        : 6.0.1
    ipywidgets       : 7.5.1
    nbformat         : 5.0.7
    traitlets        : 5.0.4
  ~~~
- macOS バージョン 10.15.7
  - Python 3.8.6
  - numpy 1.19.2
  - scipy 1.5.2
  - matplotlib 3.3.2
  - seaborn 0.11.0
  - `jupyter --version`
  ~~~
    jupyter core     : 4.6.3
    jupyter-notebook : 6.1.4
    qtconsole        : 4.7.7
    ipython          : 7.18.1
    ipykernel        : 5.3.4
    jupyter client   : 6.1.7
    jupyter lab      : not installed
    nbconvert        : 6.0.6
    ipywidgets       : 7.5.1
    nbformat         : 5.0.7
    traitlets        : 5.0.4
  ~~~

計算環境の用意が難しい方は以下のBinderからお試しいただけます（自分の計算環境での実行をお勧めします．周期解の追跡など計算にかなりの時間を要します）：

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/tmiyaji/sgc164/HEAD)
