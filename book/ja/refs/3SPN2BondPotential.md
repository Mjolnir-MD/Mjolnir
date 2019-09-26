# 3SPN2BondPotential

3SPN2ポテンシャルで使われる結合長ポテンシャルです。

{% math %}
U(v) = k(v-v_0)^2 + 100k(v - v_0)^4
{% endmath %}

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "3SPN2Bond"
topology    = "bond"
parameters  = [
    {indices = [0, 1], v0 = 1.0, k = 0.1},
    # ...
]
```

## 入力

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `v0`: 浮動小数点数型
  - 最安定点を指定します。
- `indices`: 整数の配列型
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
  - 現在`BondLength`しかこのポテンシャルをサポートしていないので、長さは2です。
