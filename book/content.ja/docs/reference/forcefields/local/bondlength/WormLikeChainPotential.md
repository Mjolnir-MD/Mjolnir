+++
title = "WormLikeChain"
weight = 600000
+++

# WormLikeChain

Worm-Like chain modelに基づいたポテンシャルです。

{{<katex display>}}
U(r) = \frac{k_B T}{p}  \left(\frac{l_c}{4} \frac{1}{1 - \frac{r}{l_c}} - \frac{r}{4} + \frac{r^2}{2l_c}\right)
{{</katex>}}

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "WormLikeChain"
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100, p = 5.0, lc = 100.0},
    # ...
]
```

## 入力

- `p`: 浮動小数点数型
  - ポリマーの持続長です。
- `lc`: 浮動小数点数型
  - ポリマーの最大長です。
- `indices`: 整数の配列型（長さ: 2）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
