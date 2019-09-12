# HardCoreExcludedVolume

硬い核と、排除体積ポテンシャルに相当する殻で構成されたポテンシャルです。

{% math %}
U(r) = \epsilon\left(\frac{\sigma}{r - r_0})^{12}
{% endmath %}

上記のパラメータ$$\sigma$$は各粒子の殻の厚さの和として定義されます。
パラメータ$$\epsilon$$は粒子によらず一定です。

## 例

```toml
[[forcefields.global]]
cutoff     = 2.5
epsilon    = 0.2
parameters = [
    {index = 0, core_radius = 3.0, soft_shell_thickness = 2.0},
    {index = 1, core_radius = 3.0, soft_shell_thickness = 2.0},
    # 必要に応じて続きます...
]
```

## 入力

- `cutoff`: 浮動小数点数型
  - カットオフ距離を決めます。
  - 相互作用している粒子ペアの殻の厚さの和に対しての相対値です。書かなければ、デフォルトの値になります。
- `epsilon`: 浮動小数点数型
  - エネルギーに相当するパラメータです。粒子に対して一律です。
- `index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。番号は0から始まります。
- `core_radius`: 浮動小数点数型
  - 柔らかい殻の厚さを表すパラメータです。冒頭の式の$$\sigma$$に相当します。
- `soft_shell_thickness`: 浮動小数点数型
  - 硬い核の半径を表すパラメータです。冒頭の式の$$r_0$$に相当します。
