+++
title = "Stoichiometric"
weight = 20000
+++

# StoichiometricInteraction

`StoichiometricInteraction`相互作用は、結合を形成する粒子同士が化学量論を保存するようにかかる相互作用です。`parameter_kinds`で指定された粒子タイプのヘテロな組み合わせのみについて適用されます。また、相互作用の化学量論係数を指定することができます。

{{<katex display>}}
V(r_{ij}) = - \epsilon \sum_{i \in a, j \in b} \frac{U(r_{ij})}{\max({\rm CO}_i, {\rm CO}_j, 1)}
{{</katex>}}

上記の式の{{<katex>}}{\rm CO}_i{{</katex>}}と{{<katex>}}{\rm CO}_j{{</katex>}}はそれぞれ

{{<katex display>}}
{\rm CO}_i = \frac{\sum_{j \in b} u(r_{ij})}{m}, {\rm CO}_j = \frac{\sum_{i \in a} u(r_{ij})}{n}
{{</katex>}}

で定義されます。ここで{{<katex>}}m{{</katex>}}は粒子bの化学量論係数、{{<katex>}}n{{</katex>}}は粒子aの化学量論係数を表します。

## 例

```toml
[[forcefields.global]]
interaction              = "Stoichiometric"
potential                = "StoichiometricUniformCubicPan"
ignore.molecule          = "Nothing"
spatial_partition.type   = {type = "CellList", margin = 0.2}
spatial_partition.margin = 0.2
particle_kinds = [
    {name = "A", coef = 1},
    {name = "B", coef = 1}
]
epsilon = 5.0  # parameter for StoichiometricInteraction
v0      = 10.0 # parameter for StoichiometricUniformCubicPanPotential
range   = 10.0 # parameter for StoichiometricUniformCubicPanPotential
paramters = [
    {index = 0, kind = "A"},
    {index = 1, kind = "A"},
    {index = 2, kind = "B"},
    {index = 3, kind = "B"},
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - `"Stoichiometric"`を指定します。
- `potential`: 文字列型
  - `U(r)`に入るポテンシャルを指定します。下記のポテンシャルが使用できます。
    - [`"StoichiometricUniformCubicPan"`]({{<relref "StoichiometricUniformCubicPanPotential.md">}})
- `ignore`: テーブルの配列型
  - この相互作用で無視する粒子ペアの条件を決めます。
  - 詳細は[GlobalForceFieldのspatial_partitionセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してください。
- `spatial_partition`: テーブル型
  - 近接リストを構築するアルゴリズムを指定します。
  - 詳細は[GlobalForceFieldのspatial_partitionセクション]({{<relref "/docs/reference/forcefields/global#ignore">}})を参照してくだい。
- `particle_kinds`: テーブルの配列型
  - この相互作用で用いる粒子の名前と対応する化学両論係数を指定します。ここで2種類の粒子のパラメータを決める必要があります。
    - `name`: 文字列型
      - `parameters`セクションで使う粒子の名前を決めます。
    - `coef`: 整数型
      - 化学両論係数を指定します。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 適用する粒子の番号を指定します。
  - `offset`: 整数型(optional、デフォルト 0.)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - `name`: 文字列型
    - `particle_kinds`で定義した粒子名を指定します。それ以外の文字列を指定した場合はエラーになります。
