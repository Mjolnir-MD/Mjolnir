# 3SPN2ExcludedVolumePotential

3SPN2系統の力場で用いられる排除体積相互作用です。

{% math %}
U(r) =
\begin{cases}
\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - 2\left(\frac{\sigma_{ij}}{r}\right)^6\right] + \epsilon & (r < \sigma_{ij})\\
0 & (r \geq \sigma_{ij})\\
\end{cases}
{% endmath %}

カットオフ距離は、常に$$\sigma$$になります。

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
    # 必要に応じて続きます...
]
```

## 入力
- `index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。番号は0から始まります。
- `kind`: 文字列型
  - 粒子の種類(`S`(sugar), `P`(phosphate), `A`, `T`, `C`, `G`)を指定します。
