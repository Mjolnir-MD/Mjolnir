+++
title = "InversePower"
weight = 500000
+++

# InversePowerPotential

Excluded Volumeポテンシャルの一般的な形です。

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{n}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "InversePower"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
ignore.groups.intra = ["chain-A"]
ignore.groups.inter = [["chain-B", "chain-C"]]

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
n          = 5
parameters = [
    {index = 0, radius = 2.0}
]
```

## 入力

- `cutoff`: 浮動小数点数型（省略可能。デフォルトで、{{<katex>}}2^{\frac{12}{n}}{{</katex>}}）
  - カットオフ長です。{{<katex>}} \sigma_{ij} {{</katex>}}との相対です。
- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを決めます。
- `n`: 整数型
  - ポテンシャルの傾きを決めます。
- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `radius`: 浮動小数点数型
  - 粒子のサイズを決めます。
  - 粒子ペアでの{{<katex>}} \sigma {{</katex>}}は`radius`の和になります。

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
