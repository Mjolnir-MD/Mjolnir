# ConstraintForceField

ConstraintForceFieldでは、特定の粒子間の距離に拘束をかけられます。

## 例

```toml
[[forcefields.constraint]]
topology      = "bond"
max_iteration = 500
tolerance     = 1e-6
parameters  = [
    {indices = [0, 1], offset = 100, v0 = 3.8},
]
```

## 入力

`ConstraintForceField`を用いる場合、`[[forcefields.constraint]]`テーブルには以下の値を設定します。

- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `max_iteration`: 整数型
  - 収束するまで最大これだけの回数補正を繰り返します。
  - 収束しなかった場合、warningが出力されますが、シミュレーションの実行は継続します。ただし、発散することが多いです。
- `tolerance`: 浮動小数点数型
  - 収束判定に用いる許容誤差を指定します。距離の差がこの値以下になれば補正の必要はないと判断されます。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - `v0`: 浮動小数点数型
    - 拘束する距離です。
