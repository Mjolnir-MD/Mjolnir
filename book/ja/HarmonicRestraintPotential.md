# HarmonicRestraint

空間中の点との距離に対して適用される調和振動子ポテンシャルです。

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

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点数型
  - 自然長を指定します。
