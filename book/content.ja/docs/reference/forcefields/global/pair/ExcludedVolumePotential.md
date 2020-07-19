+++
title = "ExcludedVolume"
weight = 400000
+++

# ExcludedVolume

排除体積効果のポテンシャルです。

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{12}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
parameters = [
    {index = 0, radius = 2.0}
]
```

## 入力

- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `index`: 整数型
  - 粒子のインデックスを指定します。
- `offset`: 整数型（省略可能）
  - インデックスのオフセットを指定します。
- `radius`: 浮動小数点数型
  - 粒子のサイズを指定します。
  - 粒子ペアでの{{<katex>}} \sigma {{</katex>}}は`radius`の和になります。
- `cutoff`: 浮動小数点数型（省略可能）
  - カットオフ長です。{{<katex>}}\sigma_{ij}{{</katex>}}に相対です。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
