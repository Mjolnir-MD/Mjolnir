# LennardJonesPotential

有名な$$(6,12)$$型のLennard-Jonesポテンシャルです。

{% math %}
U(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
{% endmath %}

ここで、上記のパラメータはローレンツ・ベルテロ則によって各ペアに対して定義されます。

{% math %}
\sigma_{ij} = \frac{\sigma_i + \sigma_j}{2},\ 
\epsilon_{ij} = \sqrt{\epsilon_i\epsilon_j}
{% endmath %}

## 例

```toml
[[forcefields.global]]
# ここにはPairInteractionなどによって要求されるパラメータが入ります。
# ...
cutoff = 2.5
parameters = [
    {index = 0, sigma = 2.0, epsilon = 10.0},
    # 必要に応じて続きます...
]
```

## 入力
- `cutoff`: 浮動小数点数型
  - カットオフ距離を決めます。
  - 相互作用している粒子ペアの半径の$$\sigma_{ij}$$に対しての相対値です。書かなければ、デフォルトの値になります。
- `index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。番号は0から始まります。
- `sigma`: 浮動小数点数型
  - 粒子径に相当するパラメータです。
- `epsilon`: 浮動小数点数型
  - エネルギーに対応するパラメータです。
