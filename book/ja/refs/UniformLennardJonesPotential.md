# UniformLennardJonesPotential

有名な$$(6,12)$$型のLennard-Jonesポテンシャルです。

[LennardJonesPotential](LennardJonesPotential.md)との唯一の違いは、パラメータが粒子に依らず一定であることです。

{% math %}
U(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
{% endmath %}

## 例

```toml
[[forcefields.global]]
# ここにはPairInteractionなどによって要求されるパラメータが入ります。
# ...
sigma      = 2.0
epsilon    = 1.0
parameters = [ # 省略可能です。省略すると、全ての粒子に適用されます。
    {index = 0},
    {index = 1},
    # 必要に応じて続きます...
]
```

## 入力

- `sigma`: 浮動小数点数型
  - 粒子径に相当するパラメータです。
- `epsilon`: 浮動小数点数型
  - エネルギーに対応するパラメータです。
- `parameters.index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。番号は0から始まります。
