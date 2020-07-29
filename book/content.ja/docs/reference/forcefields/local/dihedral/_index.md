+++
title = "DihedralAngle"
weight = 30000
bookCollapseSection = true
+++

# DihedralAngleInteraction

`DihedralAngle`相互作用は、二面角にかかる相互作用です。

以下のポテンシャルが利用可能です。

- [Cosine]({{<relref "CosinePotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [ClementiDihedral]({{<relref "ClementiDihedralPotential.md">}})
- [FlexibleLocalDihedral]({{<relref "FlexibleLocalDihedralPotential.md">}})

## 例

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Cosine"
topology    = "none"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1, 2, 3], ... },
]
```

- `interaction`: String
  - Name of the interaction. Here, it is `"BondLength"`.
- `potential`: String
  - The following potentials are available.
  - [Cosine]({{<relref "CosinePotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [ClementiDihedral]({{<relref "ClementiDihedralPotential.md">}})
  - [FlexibleLocalDihedral]({{<relref "FlexibleLocalDihedralPotential.md">}})
- `topology`: String
  - Name of the connection in [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}}).
- `parameters`: Array of Tables
  - `indices`: Array of Integers (length = 3)
    - Indices of particles that interact with each other. The index is 0-based.
  - `offset`: Integer(Optional. By default, 0.)
    - Offset of index.
  - The other parameter depends on the specified potential.

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"DihedralAngle"`を指定します。
- `potential`: 文字列型
  - [Cosine]({{<relref "CosinePotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [ClementiDihedral]({{<relref "ClementiDihedralPotential.md">}})
  - [FlexibleLocalDihedral]({{<relref "FlexibleLocalDihedralPotential.md">}})
- `topology`: 文字列型
  - [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(省略可能)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。

## Combination

いくつかの力場では、同じ二面角に複数のポテンシャルが適用されていることがあります。
そのような場合に計算を高速化するため、以下の組み合わせがサポートされています。

- `"Gaussian+FlexibleLocalDihedral"`
- `"Gaussian+Cosine"`

それぞれのポテンシャルのパラメータは、`PotentialName = {...}`のようなインラインテーブルとして設定します。

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Gaussian+FlexibleLocalDihedral"
topology    = "none"
env.ALA-ALA = [2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375]
parameters  = [
{indices = [0,1,2,3], Gaussian = {v0=-2.2, k=-0.43,sigma=0.15}, FlexibleLocalDihedral = {k=1.0, coef="ALA-ALA"}},
# ...
]
```
