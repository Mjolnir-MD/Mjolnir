# ExternalDistance

空間中の形状との距離に応じて適用される相互作用です。

## 例

```toml
[[forcefields.external]]
interaction    = "Distance"
potential      = "LennardJonesWall"
shape.name     = "AxisAlignedPlane"
shape.axis     = "+X"
shape.position = 0.0
shape.margin   = 0.5
parameters     = [
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - `"Distance"`です。
- `potential`: 文字列型
  - [`"LennardJonesWall"`](LennardJonesWallPotential.md): L-J(6,12)型のポテンシャルを用います。
  - [`"ExcludedVolumeWall"`](ExcludedVolumeWallPotential.md): 逆12乗型の斥力ポテンシャルを用います。
  - [`"ImplicitMembrane"`](ImplicitMembranePotential.md): 平面付近のエネルギーを下げるポテンシャルを用います。
- `shape`: テーブル型
  - 相互作用させる形状を指定します。

### `shape`

- `name`: 文字列型
  - 相互作用する形状の種類を決めます。
  - `"AxisAlignedPlane"`: ある軸に垂直な平面です。以下の追加パラメータを要求します。
    - `axis`: 文字列型
      - どの軸に垂直にするかを決めます。`"+X"`や`"-Z"`のように指定します。
    - `position`: 浮動小数点数型
      - 垂直な軸の上での平面の位置を指定します。
- `margin`: 浮動小数点数型
  - 相互作用する粒子のリストを管理する際のマージンの、カットオフに対する相対長さを決めます。
