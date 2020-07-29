+++
title = "RectangularBox"
weight = 20000
bookCollapseSection = true
+++

# RectangularBox

`RectangularBox`は粒子を箱の中に拘束する相互作用です。

{{<hint warning>}}
初期位置で粒子が箱の外に出ていた場合、エラーで終了します。
{{</hint>}}

{{<hint warning>}}
周期境界条件が課せられていた場合、機能しません。境界条件は設定しないようにしてください。
{{</hint>}}


## 例

```toml
[[forcefields.external]]
interaction = "RectangularBox"
potential   = "ExcludedVolumeWall"

box.lower   = [  0.0,   0.0,   0.0]
box.upper   = [100.0, 100.0, 100.0]
box.margin  = 0.4

# potential related
epsilon     = 0.1
parameters  = [
    {index = 0, radius = 1.0}, # required parameters depend on potential.
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。ここでは、`"RectangularBox"`です。
- `potential`: 文字列型
  - [`"LennardJonesWall"`]({{<relref "LennardJonesWallPotential.md">}})
  - [`"ExcludedVolumeWall"`]({{<relref "ExcludedVolumeWallPotential.md">}})
- `box`: テーブル型
  - `box.lower`: 浮動小数点数の配列型（長さ: 3）
    - 箱の座標が小さい方の頂点です。
  - `box.upper`: 浮動小数点数の配列型（長さ: 3）
    - 箱の座標が大きい方の頂点です。
  - `box.margin`: 浮動小数点数型
    - 内部で用いる近接リストのマージンです。カットオフ長に対する相対値です。
