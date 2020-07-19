+++
title = "UniformLennardJones"
weight = 200000
+++

# UniformLennardJonesPotential

よく知られたLennard-Jonesポテンシャルです。

{{<katex display>}}
U(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
{{</katex>}}

全ての粒子が同じパラメータを持っている場合のためのポテンシャルです。

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "UniformLennardJones"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}

cutoff  = 2.5
sigma   = 2.0
epsilon = 0.5
parameters = [
    {index = 0, offset = 100}, # to control which particle participates
]
```

## 入力

- `index`: 整数型
  - 粒子の番号です。最初の粒子は0番目です。
- `offset`: 整数型（省略可能）
  - 番号のオフセットです。省略可能です。グループ内番号などを使う際に便利です。
- `sigma`: 浮動小数点数型
  - 粒子のサイズを指定します。
- `epsilon`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `cutoff`: 浮動小数点数型（省略可能）
  - カットオフ長です。{{<katex>}}\sigma_{ij}{{</katex>}}に相対です。

入力の他の部分は、[Pair]({{<relref "/docs/reference/forcefields/global/pair">}})を参照してください。
