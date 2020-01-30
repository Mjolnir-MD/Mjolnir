# InversePowerPotential

距離の逆数の指定した累乗に比例する排除体積ポテンシャルです。累乗には整数のみが指定できます。
12乗を指定した場合、[ExcludedVolume](ExcludedVolumePotential.md)に相当するポテンシャルになります。

{% math %}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{n}
{% endmath %}

ここで、上記のパラメータ$$\sigma$$は各粒子の半径の和として定義されます。
パラメータ$$\epsilon$$と$$n$$は粒子によらず一定です。

## 例

```toml
[[forcefields.global]]
cutoff     = 2.5
epsilon    = 0.2
n          = 5
parameters = [
    {index = 0, radius = 2.0}
]
```


## 入力

- `cutoff`: 浮動小数点数型
  - カットオフ距離を決めます。
  - 相互作用している粒子ペアの半径の和に対しての相対値です。書かなければ、デフォルトの値になります。デフォルトの値は$$2^{\frac{12}{n}}$$です。
- `epsilon`: 浮動小数点数型
  - エネルギーに相当するパラメータです。粒子に対して一律です。
- `n`: 整数型
  - 粒子の硬さに相当するパラメータです。粒子に対して一律です。
- `radius`: 浮動小数点数型
  - 粒子系に相当するパラメータです。
