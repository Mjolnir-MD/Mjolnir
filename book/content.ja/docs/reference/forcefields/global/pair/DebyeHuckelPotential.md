+++
title = "DebyeHückel"
weight = 300000
+++

# DebyeHückelPotential

Debye-Hückelの式に基づいた、イオン溶液内での静電相互作用のモデルです。

{{<katex display>}}
U(r_{ij}) = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_k r_{ij}} \exp(-r_{ij}/\lambda_D)
{{</katex>}}

{{<katex display>}}
\lambda_D = \sqrt{\frac{\epsilon_0\epsilon_k}{2\beta N_A e_c^2 I}}
{{</katex>}}

温度・イオン強度に依存した水の電気伝導度は以下の論文の式に従います。

- Sambriski, E. J. et al., (2009) Biophys. J.
- Hinckley, D. M. et al., (2013) JCP.
- Freeman, G. S., et al., (2014) JCP.

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "DebyeHuckel"
ignore.particles_within.bond = 3
spatial_partition.type = "CellList"
spatial_partition.margin = 0.2
cutoff     = 5.5
parameters = [
    {index = 0, charge = 1.0},
]
```

## 入力

このポテンシャルは、[system attribute]({{<relref "/docs/reference/system#attribute">}})で`temperature`と `ionic_strength`を設定していることを要求します。

- `cutoff`: 浮動小数点数型（省略可能。デフォルトでは5.5）
  - カットオフ長です。デバイ長との相対です。
- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `charge`: 浮動小数点数型
  - 粒子の電荷です。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
