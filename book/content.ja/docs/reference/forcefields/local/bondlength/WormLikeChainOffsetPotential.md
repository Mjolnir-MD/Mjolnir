+++
title = "WormLikeChainOffset"
weight = 900000
+++

# WormLikeChainOffset

Worm-Like chain modelに基づいた、二粒子間の距離に関するオフセットを考慮したポテンシャルです。

{{<katex display>}}
U(r) = \begin{cases}
0 & (r < l_0) \\
 \frac{k_B T}{p}  \left(\frac{l_c}{4} \left[ \frac{1}{1 - \frac{r - l_0}{l_c}} - 1 \right] - \frac{r - l_0}{4} + \frac{\left(r - l_0\right)^2}{2l_c}\right) & (r \geq l_0) \\
\end{cases}
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "WormLikeChainOffset"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, p = 5.0, lc = 100.0, l0 = 30.0},
    # ...
]
```

## 入力

- `p`: 浮動小数点数型
  - ポリマーの持続長です。
- `lc`: 浮動小数点数型
  - ポリマーの最大長です。
- `l0`: 浮動小数点数型
  - 二粒子間の距離に対してかかるオフセットです。
- `indices`: 整数の配列型（長さ: 2）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
