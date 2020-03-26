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
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
integrator.parameters = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ... 必要に応じて続きます
]
```

## 入力

`delta_t`などの他のパラメータは[Simulator](Simulator.md)で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"UnderdampedLangevin"`です。
- `remove`: テーブル型 (optional)
  - 系全体の並進・回転速度成分を取り除くことができます。`true`にすると取り除かれます。
  - `rescale`が`true`になっていた場合、全体の速度ベクトルをリスケールすることで速度を減算した分の運動エネルギーを補填します。
  - 省略した場合、全て`false`になります。
- `parameters`: テーブルの配列型
  - 各粒子の摩擦係数$$\gamma_i$$を指定します。

{% hint style='danger' %}
バージョン1.5.xまではここで乱数シードを設定していましたが、1.6.x以降から各
`[simulator]`で設定するようになりました。`simulator.integrator.seed`は
1.7.x以降削除予定です。
{% endhint %}
