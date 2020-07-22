# WormLikeChainPotential

WormLikeChainモデルを元にしたポテンシャルです。持続長と最大伸長時の長さからポリマーの伸長性を表現します。

{% math %}
U(l) = \frac{k_B T}{p} \left( \frac{lc}{4} \left[ \frac{1}{1 - \frac{l}{lc}} - 1 \right] - \frac{l}{4} + \frac{l^2}{2lc} \right)
{% endmath %}

## 例

```toml
[[forcefields.local]]
# ここにはBondLengthなどによって要求されるパラメータが入ります。
# ...
parameters = [
    {indices = [0, 1], p = 3.8, lc = 390.0},
    # 必要に応じて続きます...
]
```

## 入力

- `p`: 浮動小数点数型
  - ポリマーの持続長に相当するパラメータを指定します。
- `lc`: 浮動小数点数型
  - ポリマーの最大伸長時の長さに相当するパラメータを指定します。
- `indices`: 整数の配列型
  - どの粒子の間に適用するかを指定します。最初の粒子は0番目です。
