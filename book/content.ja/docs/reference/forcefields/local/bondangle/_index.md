+++
title = "BondAngle"
weight = 20000
bookCollapseSection = true
+++

# BondAngleInteraction

`BondAngle`相互作用は、結合角にかかる相互作用です。

以下のポテンシャルが利用可能です。

- [Harmonic]({{<relref "HarmonicPotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [FlexibleLocalAngle]({{<relref "FlexibleLocalAnglePotential.md">}})

## 例

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Harmonic"
topology    = "none"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1, 2], offset = 100, ... },
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"BondAngle"`を指定します。
- `potential`: 文字列型
  - 以下のポテンシャルが利用可能です。
  - [Harmonic]({{<relref "HarmonicPotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [FlexibleLocalAngle]({{<relref "FlexibleLocalAnglePotential.md">}})
- `topology`: 文字列型
  - [`"Topology"`]({{<relref "/docs/reference/forcefields/Topology.md">}})での名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型（長さ: 3）
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。
