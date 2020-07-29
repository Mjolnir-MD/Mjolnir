+++
title = "BondLength"
weight = 10000
bookCollapseSection = true
+++

# BondLengthInteraction

`BondLength`相互作用は、その名の通り結合長にかかる相互作用です。

以下のポテンシャルが利用可能です。

- [Harmonic]({{<relref "HarmonicPotential.md">}})
- [Gaussian]({{<relref "GaussianPotential.md">}})
- [GoContact]({{<relref "GoContactPotential.md">}})
- [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
- [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
- [WormLikeChain]({{<relref "WormLikeChainPotential.md">}})
- [3SPN2Bond]({{<relref "3SPN2BondPotential.md">}})

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1], offset = 100, ... },
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"BondLength"`を指定します。
- `potential`: 文字列型
  - 利用可能なポテンシャルは以下のとおりです。
  - [Harmonic]({{<relref "HarmonicPotential.md">}})
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [GoContact]({{<relref "GoContactPotential.md">}})
  - [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
  - [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
  - [WormLikeChain]({{<relref "WormLikeChainPotential.md">}})
  - [3SPN2Bond]({{<relref "3SPN2BondPotential.md">}})
- `topology`: 文字列型
  - [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})における結合の名称です。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型(長さ: 2)
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(省略可能)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。
