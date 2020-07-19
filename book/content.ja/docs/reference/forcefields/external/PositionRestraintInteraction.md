+++
title = "PositionRestraint"
weight = 10000
+++

# PositionRestraint

粒子を空間中の点に、あるいは空間中のある点から一定距離に保つ相互作用です。

## 例

```toml
[[forcefields.external]]
interaction = "PositionRestraint"
potential   = "Harmonic"
parameters  = [
    {index = 0, position = [0.0, 0.0, 0.0], k = 0.1, v0 = 10.0},
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"PositionRestraint"`を指定します。
- `potential`: 文字列型
  - 点との距離に適用するポテンシャルを指定します。
  - [`"Harmonic"`](HarmonicRestraintPotential.md): 調和振動子ポテンシャルを使用します。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 適用する粒子の番号を指定します。
  - `offset`: 整数型（省略可能）
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - `position`: 浮動小数点数の配列型
    - 粒子を固定する位置を指定します。
  - `k`: 浮動小数点数型
    - ポテンシャルの強さを決めます。
  - `v0`: 浮動小数点数型
    - 点からの距離を決めます。
    - 点の位置に固定したい場合、0にしてください。
