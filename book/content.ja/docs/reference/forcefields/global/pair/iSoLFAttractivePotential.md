+++
title = "iSoLF"
weight = 100000
+++

# iSoLFAttractivePotential

iSoLFは以下の文献で開発された粗視化膜モデルで、iSoLFAttractiveはその引力ポテンシャルです。

- Diego Ugarte La Torre and Shoji Takada (2020) J. Chem. Phys 153, 205101
  https://doi.org/10.1063/5.0026342

{{<katex display>}}
U(r) =
\begin{cases}
-\epsilon_{ij} & r_{ij} < \sqrt[6]{2}\sigma_{ij}\\
-\epsilon_{ij} \cos^2\left[\frac{\pi}{2\omega_{ij}}(r_{ij} - \sqrt[6]{2}\sigma_{ij}) \right] & (\sqrt[6]{2}\sigma_{ij} < r_{ij} < \sqrt[6]{2}\sigma_{ij} + \omega_{ij})\\
0 & (\sqrt[6]{2}\sigma_{ij} + \omega_{ij} < r_{ij})
\end{cases}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "iSoLFAttractive"
ignore.particles_within = {bond = 1, angle = 1}
ignore.group.inter = [
    ["T1", "T3"]
]
spatial_partition = {type = "CellList", margin = 0.5}
env.popc_epsilon = 0.416
env.popc_omega   = 9.867
env.popc_sigma_T = 7.111
parameters = [
{index =     2, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     3, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     4, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     7, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     8, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
{index =     9, sigma = "popc_sigma_T", epsilon = "popc_epsilon", omega = "popc_omega"},
# ...
]
```

## 入力

{{<katex>}} \sigma {{</katex>}}と{{<katex>}} \epsilon {{</katex>}}と{{<katex>}} \omega {{</katex>}}の計算には、Lorentz-Berthelot則が用いられます。

{{<katex display>}}
\sigma_{ij}   = \frac{\sigma_i + \sigma_j}{2} \\
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j} \\
\omega_{ij}   = \frac{\omega_i\omega_j}{2}
{{</katex>}}

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `sigma`: 浮動小数点数型
  - 粒子のサイズを指定します。
- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `omega`: 浮動小数点数型
  - ポテンシャルの届く範囲を指定します。

このポテンシャルは{{<katex>}} r = \sqrt[6]{2}\sigma_{ij} + \omega_{ij} {{</katex>}}で正確に0となるため、カットオフには常にこの値が使用されます。
よって、`cutoff`を指定する必要はありません。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
