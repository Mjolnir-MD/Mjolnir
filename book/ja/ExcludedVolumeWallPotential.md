# ExcludedVolumeWall

空間中の形状と粒子の間に適用される、逆12乗の斥力ポテンシャルです。

{% math %}
U(x) = \epsilon\left(\frac{r_0}{x}\right)^{12}
{% endmath %}

## 例

`shape`パラメータについては[ExternalDistanceInteraction](ExternalDistanceInteraction.md)を参照して下さい。

```toml
[[forcefields.external]]
interaction    = "Distance"
shape.name     = "AxisAlignedPlane"
shape.axis     = "+X"
shape.position = 0.0
shape.margin   = 0.5

potential   = "ExcludedVolume"
epsilon     = 0.1
parameters  = [
    {index = 0, radius = 1.0},
    # ...
]
```

## 入力

- `epsilon`: 浮動小数点数型
  - エネルギーに相当するパラメータです。粒子によらず一定です。
- `radius`: 浮動小数点数型
  - 粒子径に相当するパラメータです。

