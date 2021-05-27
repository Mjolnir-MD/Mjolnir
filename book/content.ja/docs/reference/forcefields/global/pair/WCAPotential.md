+++
title = "WCA"
weight = 700000
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

パラメータを指定するには二通りの方法があります。

個々の粒子のパラメータを指定し、Lorentz-Berthelot則でペアのパラメータを計算する方法と、

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

粒子の名前を指定し、表でペアのパラメータを与える方法があります。

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJonesAttractive"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

table.A.A = {sigma = 1.0, epsilon = 2.0}
table.A.B = {sigma = 3.0, epsilon = 1.0} # B.A will be the same
table.B.B = {sigma = 1.0, epsilon = 1.5}
parameters = [
    {index = 0, offset = 100, name = "A"},
    {index = 1, offset = 100, name = "B"},
    # ...
]
```

## 入力

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `sigma`: 浮動小数点数型
  - 粒子のサイズを指定します。テーブルを用意する場合、これは不要です。
- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。テーブルを用意する場合、これは不要です。
- `name`: 文字列型
  - 粒子の名前です。テーブルを用意する場合、これは必須です。

このポテンシャルは{{<katex>}} r = \sigma\sqrt[6]{2} {{</katex>}}で正確に0となるため、カットオフには常にこの値が使用されます。
よって、`cutoff`を指定する必要はありません。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
