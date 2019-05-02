# ClementiDihedral

Clementi-Goポテンシャルで用いられる二面角ポテンシャルです。

{% math %}
U(v) = k_1(1-\cos(v-v_0)) + k_3(1-\cos(3(v-v_0)))
{% endmath %}

以下の論文で提案された力場です。

- C. Clementi, H. Nymeyer, J. Onuchic, (2000) JMB

## 例

```toml
[[forcefields.local]]
parameters = [
    {indices = [0,1,2,3], v0 = -2.20, k1 = 1.0, k3 = 0.5},
    # ...
]
```

## 入力

- `v0`: 浮動小数点数型
  - 最安定角度を指定します。
- `k1`: 浮動小数点数型
  - ポテンシャルの強さを決めるパラメータの一つです。
- `k3`: 浮動小数点数型 
  - ポテンシャルの強さを決めるパラメータの一つです。
- `indices`: 整数の配列型
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
  - `ClementiDihedral`は二面角相互作用のためのポテンシャルなので、4つ指定します。
