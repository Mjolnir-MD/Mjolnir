+++
title = "3SPN2Base"
weight = 20000
+++

# 3SPN2BaseBaseInteraction

`3SPN2BaseBaseInteraction`は3SPN2系粗視化DNA力場の相互作用です。

以下のポテンシャルが利用可能です。
- `"3SPN2"`: Hinckley et al., (2013) JCP
- `"3SPN2C"`: Freeman et al., (2014) JCP

## 例

```toml
[[forcefields.global]]
interaction = "3SPN2BaseBase"
potential   = "3SPN2"
ignore.particles_within.nucleotide = 3
spatial_partition = {type = "CellList", margin = 0.2}
parameters  = [
# `nucleotide` index starts from 5' and ends at 3'.
{strand = 0, nucleotide =  0,          S =   0, B =   1, offset = 100, Base = "A"},
{strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, offset = 100, Base = "T"},
{strand = 0, nucleotide =  2, P =   5, S =   6, B =   7, offset = 100, Base = "C"},
# ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の種類を設定します。`3SPN2BaseBase`です。
- `potential`: 文字列型
  - 以下のポテンシャルが利用可能です。
  - `3SPN2` : 3SPN.2ポテンシャルのパラメータを使います。
  - `3SPN2C`: 3SPN.2Cポテンシャルのパラメータを使います。
- `topology`: 文字列型
  - [`"Topology"`]({{<relref "/docs/reference/forcefields/Topology.md">}})に設定する名前を指定します。
  - 3SPN2力場では3つのヌクレオチドが間にないと塩基対を形成しないので、`BaseStacking`の
    トポロジーとして設定した結合(ここでは`"nucleotide"`)を3つ分無視して下さい。
  - 詳細は[GlobalForceFieldのignoreセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してください。
- `spatial_partition`: Table
  - 近接リストを構築するアルゴリズムを指定します。
  - 詳細は[GlobalForceFieldのspatial_partitionセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してください。
- `parameters`:テーブルの配列型
  - `strand`: 整数型
    - ストランドの番号です。0-basedなインデックスです。
  - `nucleotide`: Integer
    - ヌクレオチドの番号です。0-basedなインデックスです。
  - `P, S, B`: 整数型
    - リン酸(P), 糖(S), 塩基(B)の番号です。0-basedなインデックスです。
  - `Base`: String
    - `"A"`, `"T"`, `"C"`, `"G"`のいずれかです。
  - `offset`: 整数型(optional)
    - 各粒子の粒子番号に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 端のヌクレオチドには、多くの場合リン酸がありません。
