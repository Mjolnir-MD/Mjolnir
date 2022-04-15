+++
title = "Pair"
weight = 10000
bookCollapseSection = true
+++

# GlobalPair

`GlobalPairInteraction`はパラメータを持つ全ての粒子の間に働く相互作用です。

以下のポテンシャルが利用可能です。

- [LennardJones]({{<relref "LennardJonesPotential.md">}})
- [UniformLennardJones]({{<relref "UniformLennardJonesPotential.md">}})
- [DebyeHuckel]({{<relref "DebyeHuckelPotential.md">}})
- [ExcludedVolume]({{<relref "ExcludedVolumePotential.md">}})
- [InversePower]({{<relref "InversePowerPotential.md">}})
- [HardCoreExcludedVolume]({{<relref "HardCoreExcludedVolumePotential.md">}})
- [WCA]({{<relref "WCAPotential.md">}})
- [iSoLFAttractive]({{<relref "iSoLFAttractivePotential.md">}})
- [3SPN2ExcludedVolume]({{<relref "3SPN2ExcludedVolumePotential.md">}})
- [UniformCubicPan]({{<relref "UniformCubicPanPotential.md">}})

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
ignore.groups.intra = ["chain-A"]
ignore.groups.inter = [["chain-B", "chain-C"]]
spatial_partition.type = "CellList"
spatial_partition.margin = 0.2
parameters = [
    {index = 0, offset = 100, ...}, # required parameter depends on the potential.
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。ここでは、`"Pair"`です。
- `potential`: 文字列型
  - 以下のポテンシャルが利用可能です。
  - [`"LennardJones"`]({{<relref "LennardJonesPotential.md">}})
  - [`"UniformLennardJones"`]({{<relref "UniformLennardJonesPotential.md">}})
  - [`"DebyeHuckel"`](DebyeHuckelPotential.md)
  - [`"ExcludedVolume"`](ExcludedVolumePotential.md)
  - [`"InversePower"`](InversePowerPotential.md)
  - [`"HardCoreExcludedVolume"`](HardCoreExcludedVolumePotential.md)
- `ignore`: テーブル型
  - この相互作用で無視する粒子ペアの条件を決めます。
  - 詳細は[GlobalForceFieldのignoreセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してください。
- `spatial_partition`: Table
  - 近接リストを構築するアルゴリズムを指定します。
  - 詳細は[GlobalForceFieldのspatial_partitionセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してください。
- `parameters`: Array of Tables
  - `index`: 整数型
    - 粒子の番号です。最初の粒子は0番目です。
  - `offset`: 整数型 (省略可能)
    - 各粒子の粒子番号に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
