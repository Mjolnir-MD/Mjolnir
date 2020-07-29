+++
title = "HardCore"
weight = 600000
+++

# HardCoreExcludedVolume

中心に重なり合わない硬い核を持つ排除体積ポテンシャルです。

{{<katex display>}}
U(r) = \epsilon\left(\frac{\sigma}{r - r_0}\right)^{12}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "HardCoreExcludedVolume"

ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1

spatial_partition.type = "CellList"
spatial_partition.margin = 0.2

cutoff     = 2.5
epsilon    = 0.2
parameters = [
    {index = 0, core_radius = 3.0, soft_shell_thickness = 2.0},
    {index = 1, core_radius = 3.0, soft_shell_thickness = 2.0},
]
```

## 入力

- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `index`: 整数型
  - 粒子のインデックスを指定します。
- `offset`: 整数型（省略可能）
  - インデックスのオフセットを指定します。
- `core_radius`: 浮動小数点数型
  - 核の半径です。上の式での{{<katex>}} r_0 {{</katex>}}に相当します。
  - 粒子ペアでの{{<katex>}} r_0 {{</katex>}}は`core_radius`の和になります。
- `soft_shell_thickness`: 浮動小数点数型
  - 核の外側の半径です。上の式での{{<katex>}} \sigma {{</katex>}}に相当します。
  - 粒子ペアでの{{<katex>}}\sigma_{ij}{{</katex>}}は`soft_shell_thickness`の和になります。
- `cutoff`: 浮動小数点数型（省略可能）
  - カットオフ長です。{{<katex>}}\sigma_{ij}{{</katex>}}に相対です。

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
