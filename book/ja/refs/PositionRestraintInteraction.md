# PositionRestraint

空間中の点との距離に応じて適用される相互作用です。

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
  - `"PositionRestraint"`を指定します。
- `potential`: 文字列型
  - 点との距離に適用するポテンシャルを指定します。
  - [`"Harmonic"`](HarmonicRestraintPotential.md): 調和振動子ポテンシャルを使用します。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 適用する粒子の番号を指定します。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - `position`: 浮動小数点数の配列型
    - 粒子を固定する位置を指定します。
- その他のパラメータはポテンシャルに依存します。
