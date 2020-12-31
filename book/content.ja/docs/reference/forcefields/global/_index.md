+++
title = "Global"
weight = 2000
bookCollapseSection = true
+++

# GlobalForceField

`GlobalForceField`は相互作用に参加する粒子全ての間にかかる力です。

## [GlobalPair]({{<relref "pair">}})

2粒子間の距離に依存する相互作用です。

## [3SPN2BaseBase](3SPN2BaseBaseInteraction.md)

3SPN2粗視化DNAモデルで用いられる相互作用です。

## [ProteinDNANonSpecific](ProteinDNANonSpecificInteraction.md)

水素結合の粗視化モデルです。

## 相互作用間で共通した部分

全ての`[[forcefields.global]]`テーブルは、相互作用ペアから除外する条件を記述する`ignore`と、近接リストを作るアルゴリズムを記述する`spatial_partition`を共通して持ちます。

### ignore

除外する条件としては、3種類の条件を指定できます。

まず`particles_within`では、[Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})上で結合を辿って見つかるペアの相互作用を除外することができます。

続いて`group`では、[System]({{<relref "/docs/reference/system">}})で各粒子に定義した`group`に基づいて、`group`間・`group`内での相互作用を一気に除外することができます。

最後に`molecule`は、`bond`を伝ってたどり着ける粒子の集合を`molecule`とし、`molecule`間・`molecule`内の相互作用を一気に除外することができます。

この設定は全て省略可能です。省略した場合、どのペアも無視せず、全ての相互作用を計算します。

- `ignore.particles_within`: テーブル型（省略可能）
  - 一定個数の結合を辿って見つかる範囲にある粒子のペアを除外します。
  - このテーブルで定義されているキーがたどる結合の名前で、値はたどる最大回数です。
  - 例えば、`ignore.particles_within.bond = 3`だと`bond`という名前を[Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})上で3回までたどって見つかる粒子のペアを除外します。
 ignores particles that are connected via 1, 2, or 3 consecutive bond connections.
- `ignore.group`: テーブル型（省略可能）
  - `group`内・間のペアを除外します。
  - `intra`: 文字列の配列型
    - グループ**内**のペアを除外します。適用するグループのリストを渡します。
  - `inter`: 文字列の配列の配列型
    - グループ**間**のペアを除外します。適用するグループのリストを渡します。
  - あとで例が示されます。
- `ignore.molecule`: 文字列型（省略可能）
  - `molecule`間・内の相互作用を除外します。
  - `"Nothing"`: 何も除外しません。全てのペアを計算します。
  - `"Self"`: 同じ`molecule`にあるペアは除外します。
  - `"Others"`: 異なる`molecule`にあるペアは除外します。

以下のような系があるとします。

```toml
[[systems]]
particles = [
    {... group = "A"}, # particle 0
    {... group = "A"}, #          1
    {... group = "A"}, #          2
     
    {... group = "B"}, # particle 3
    {... group = "B"}, #          4
    {... group = "B"}, #          5
     
    {... group = "C"}, # particle 6
    {... group = "C"}, #          7
    {... group = "C"}, #          8
]
```

以下の相互作用は、

```toml
[[forcefields.global]]
interaction = "Pair"
ignore.group.intra = ["A", "B", "C"]
```

粒子0と1、0と2、1と2は、同じグループに属しているペアなので除外されますが、0と3、0と6、3と6の間の相互作用は、異なるグループに属している相互作用なので、除外れず計算されます。

また、以下の場合は、

```toml
[[forcefields.global]]
interaction = "Pair"
ignore.group.inter = [["A", "B"], ["A", "C"]]
```

粒子0と3, 0と4, 0と5はグループAとBの間のペアなので、また3と6, 3と7, 3と8はグループAとCの間のペアなので除外されますが、0と1、3と4、3と6は計算されます。

### spatial_partition

近接リストを生成するためのアルゴリズムを指定します。

- `type`: 文字列型
  - アルゴリズムの種類を決めます。以下が使用可能です。
  - `"CellList"`: セルリストから近接リストを作ります。
  - `"RTree"`: R木を構築し、近接リストを作ります。
  - `"VerletList"`: 全ペアを見て近接リストを作ります。相互作用に関わる粒子の数が十分少ない場合は拘束なことがあります。
- `margin`: 浮動小数点数型
  - 近接リストで用いるマージンの長さです。カットオフ長との相対値です。
  - これは実行効率には影響しますが、精度には影響しません。最適な値はポテンシャルと系に依存します。
