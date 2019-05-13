# UnderdampedLangevin

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

{% math %}
m_i\boldsymbol{a}_i = \boldsymbol{f}_i - m_i\gamma_i\boldsymbol{v}_i + m_i\xi_i
{% endmath %}

以下の論文で提案されている手法です。

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

また、以下のシミュレータで推奨されている手法です。

- H. Kenzaki et al., (2011) JCTC.

## 例

```toml
[simulator]
integrator.type = "UnderdampedLangevin"
integrator.seed = 12345
integrator.parameters = [
    {index = 0, gamma = 1.0},
    {index = 1, "γ"   = 1.0},
    # ... 必要に応じて続きます
]
```

## 入力

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"UnderdampedLangevin"`です。
- `seed`: 整数型
  - 乱数のシードを設定します。
- `parameters`: テーブルの配列型
  - 各粒子の摩擦係数$$\gamma_i$$を指定します。`gamma`と`"γ"`のどちらでも構いません。
