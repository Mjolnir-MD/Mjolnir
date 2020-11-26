+++
title = "WCA"
weight = 100000
+++

# WCAPotential

よく知られたWCA (Weeks-Chandler-Andersen) ポテンシャルです。Lennard-Jonesの斥力項のみを取り出した形のポテンシャルです。

{{<katex display>}}
U(r) =
\begin{cases}
4\epsilon \left[\left(\frac{\sigma_{ij}}{r}\right)^{12} - \left(\frac{\sigma_{ij}}{r}\right)^6\right] + \epsilon & (r < \sigma_{ij}\sqrt[6]{2})\\
0 & (r \geq \sigma_{ij}\sqrt[6]{2})\\
\end{cases}
{{</katex>}}

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "WCA"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

parameters = [
    {index = 0, offset = 100, sigma = 2.0, epsilon = 10.0},
]
```

## 入力

{{<katex>}} \sigma {{</katex>}} と {{<katex>}} \epsilon {{</katex>}}の計算には、Lorentz-Berthelot則が用いられます。

{{<katex display>}}
\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2} \\
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j}
{{</katex>}}

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `sigma`: 浮動小数点数型
  - 粒子のサイズを指定します。
- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。

このポテンシャルは{{<katex>}} r = \sigma\sqrt[6]{2} {{</katex>}}で正確に0となるため、カットオフには常にこの値が使用されます。
よって、`cutoff`を指定する必要はありません。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
