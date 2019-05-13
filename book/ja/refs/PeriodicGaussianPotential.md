# PeriodicGaussianPotential

このポテンシャルは、二面角相互作用としてガウシアンポテンシャルを使うために、
最安定点との差を境界を考慮して計算するポテンシャルです。

それ以外は、[Gaussian](GaussianPotential.md)ポテンシャルと違いはありません。

{% math %}
U(v) = k\exp\left(\frac{-(v - v_0)^2}{2\sigma^2}\right)
{% endmath %}

## 例

```toml
[[forcefields.local]]
# ここにはBondLengthなどによって要求されるパラメータが入ります。
# ...
parameters = [
    {indices = [0, 1], v0 = 1.0, k = 100.0, sigma = 5.0},
    {indices = [0, 1], v0 = 1.0, k = 100.0, "σ"   = 5.0},
    # 必要に応じて続きます...
]
```

## 入力

上の例では対応する粒子の番号が2つ書かれていますが、この個数は相互作用の種類によって変動します。

- `k`: 浮動小数点数型
  - ポテンシャルの強さを指定します。
- `sigma`: 浮動小数点数型
  - ポテンシャルが効果を及ぼす幅を指定します。
  - `sigma`も、`"σ"`(GREEK SMALL LETTER SIGMA)も用いることができます。
- `v0`: 浮動小数点数型
  - 最安定点を指定します。
