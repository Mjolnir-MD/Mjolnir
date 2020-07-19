+++
title = "3SPN2ExcludedVolume"
weight = 700000
+++

# 3SPN2ExcludedVolumePotential

3SPN2粗視化DNAモデルの排除体積ポテンシャルです。

{{<katex display>}}
U(r) =
\begin{cases}
\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - 2\left(\frac{\sigma_{ij}}{r}\right)^6\right] + \epsilon & (r < \sigma_{ij})\\
0 & (r \geq \sigma_{ij})\\
\end{cases}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential = "3SPN2ExcludedVolume"
ignore.particles_within.bond       = 1
ignore.particles_within.angle      = 1
ignore.particles_within.dihedral   = 1
ignore.particles_within.nucleotide = 1
spatial_partition = {type = "CellList", margin = 0.4}
parameters = [
    {index =   0, kind = "S"},
    {index =   1, kind = "A"},
    {index =   2, kind = "P"},
    {index =   3, kind = "S"},
]
```

## 入力

- `index`: 整数型
  - 粒子のインデックスを指定します。
- `offset`: 整数型（省略可能）
  - インデックスのオフセットを指定します。
- `kind`: 文字列
  - 以下のどれかです。
  - `"S"`(sugar)
  - `"P"`(phosphate)
  - `"A"`(アデニン塩基).
  - `"T"`(チミン塩基).
  - `"C"`(シトシン塩基).
  - `"G"`(グアニン塩基).

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
