# RecutangularBox

空間中の箱の内部に粒子を留めるための相互作用です。

{% hint style='info' %}
初期状態で粒子が箱の外に出ていた場合、エラーを表示して終了します。
また、周期境界条件が適用されていた場合も、エラーを表示して終了します。
{% endhint %}

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
    {index = 0, radius = 1.0},
]
```

## 入力

- `interaction`: 文字列型
  - `"RectangularBox"`です。
- `potential`: 文字列型
  - [`"LennardJonesWall"`](LennardJonesWallPotential.md): L-J(6,12)型のポテンシャルを用います。
  - [`"ExcludedVolumeWall"`](ExcludedVolumeWallPotential.md): 逆12乗型の斥力ポテンシャルを用います。
- `box`: テーブル型
  - `box.lower`: 浮動小数点数の配列型
    - 箱の下限の座標を決めます。
  - `box.upper`: 浮動小数点数の配列型
    - 箱の上限の座標を決めます。
  - `box.margin`: 浮動小数点数型
    - 箱の上限の座標を決めます。

その他のパラメータはポテンシャルの種類に依存します。詳細は各ポテンシャルのページを参照してください。
