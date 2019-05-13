# GoContactPotential

Go-コンタクトポテンシャルは以下のような形のポテンシャルです。

{% math %}
U(r) = k\left[5\left(\frac{r_0}{r}\right)^{12} - 6\left(\frac{r_0}{r}\right)^{10}\right]
{% endmath %}

粗視化分子動力学において用いられる、粒子間の距離に対して使用されるポテンシャルです。
[`BondLength`](BondLengthInteraction.md)または[`Contact`](GoContactPotential.md)相互作用として使います。

## 例

```toml
[[forcefields.local]]
# ここにはBondLengthなどによって要求されるパラメータが入ります。
# ...
parameters = [
    {indices = [0, 1], v0 = 1.0, k = 0.1},
    # 必要に応じて続きます...
]
```

## 入力

- `v0`: 浮動小数点数型
  - 最安定距離を指定します。
- `k`: 浮動小数点数型
  - パラメータの強さを指定します。
- `indices`: 整数の配列型
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
  - `GoContact`は距離に関連した相互作用でしか使用しないので、常に2個指定します。
