+++
title = "GoContact"
weight = 300000
+++

# GoContactPotential

Go-コンタクトポテンシャルは以下のような形のポテンシャルです。

{{<katex display>}}
U(r) = k\left[5\left(\frac{r_0}{r}\right)^{12} - 6\left(\frac{r_0}{r}\right)^{10}\right]
{{</katex>}}

粗視化分子動力学の構造依拠モデルで、参照構造で距離が近かった粒子同士のコンタクトとして用いられるポテンシャルです。

## 例

```toml
[[forcefields.local]]
interaction = "Contact"
potential   = "GoContact"
topology    = "contact"
parameters = [
    {indices = [0, 1], v0 = 1.0, k = 0.1},
    # ...
]
```

## 入力

- `v0`: 浮動小数点数型
  - 最安定距離を指定します。
- `k`: 浮動小数点数型
  - パラメータの強さを指定します。
- `indices`: 整数の配列型（長さ: 2）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。
