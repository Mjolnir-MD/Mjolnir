+++
title = "Topology"
weight = 6000
+++

# Topology

`LocalForceField`は`GlobalForceField`に影響することがあります。
例えば、結合長ポテンシャルが適用されている粒子同士には排除体積効果が働かない、というケースがあります。

この情報を`LocalForceFiled`と`GlobalForceField`と共有するため、`Topology`クラスがあり、`System`に格納されています。
`Topology`クラスには粒子をノード、`LocalForceField`をエッジとしたグラフ構造が格納されており、その上を検索することができます。

`LocalForceField`では、それぞれに名前を付けて`Topology`として登録することができます。
`GlobalForceField`では、名前のついたエッジを辿って見つかる粒子を無視することができます。

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [0, 1], ... },
    # ...
]

[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"
ignore.particles_within.bond    = 3 # ignore particles within 3 bonds.
ignore.particles_within.contact = 1 # ignore particles within 1 contact.
# ...
parameters = [
    # ...
]
```
