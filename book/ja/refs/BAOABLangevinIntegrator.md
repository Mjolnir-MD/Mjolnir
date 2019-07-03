# BAOABLangevin

ランジュバン方程式に従い、温度・体積・粒子数 (NVT) 一定のシミュレーションを行います。

{% math %}
m_i\boldsymbol{a}_i = \boldsymbol{f}_i - m_i\gamma_i\boldsymbol{v}_i + m_i\xi_i
{% endmath %}

以下の論文で提案された手法です。

- Leimkuhler B, Matthews C. Appl. Math. Res. Exp. (2013)
- Leimkuhler B, Matthews C. J. Chem. Phys. (2013)

## 例

```toml
[simulator]
integrator.type = "BAOABLangevin"
integrator.seed = 12345
integrator.parameters = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ... 必要に応じて続きます
]
```

## 入力

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"BAOABLangevin"`です。
- `seed`: 整数型
  - 乱数のシードを設定します。
- `parameters`: テーブルの配列型
  - それぞれの粒子について、摩擦係数$$\gamma_i$$を指定します。