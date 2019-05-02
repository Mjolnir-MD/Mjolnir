# ExcludedVolume

距離の逆12乗に比例する排除体積ポテンシャルです。

{% math %}
U(r) = \epsilon\left(\frac{\sigma}{r}\right)^{12}
{% endmath %}

ここで、上記のパラメータ$$\sigma$$は各粒子の半径の和として定義されます。

## 例

```toml
[[forcefields.global]]
epsilon    = 0.2
parameters = [
    {index = 0, radius = 2.0}
]
```

## 入力

- `epsilon`: 浮動小数点数型
  - エネルギーに相当するパラメータです。粒子に対して一律です。
- `index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。
- `radius`: 浮動小数点数型
  - 粒子径に相当するパラメータです。

