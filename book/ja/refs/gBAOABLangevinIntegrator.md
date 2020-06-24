# gBAOABLangevin

ランジュバン方程式に従い、温度・体積・粒子数 (NVT) 一定のシミュレーションを行います。
`BAOABLangevin`と異なり、粒子同士の距離拘束を取り扱うことができます。

以下の論文で提案された手法です。

- Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)

## 例

```toml
[simulator]
integrator.type = "g-BAOABLangevin"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ... 必要に応じて続きます
]
```

## 入力

`delta_t`などの他のパラメータは[Simulator](Simulator.md)で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"g-BAOABLangevin"`です。
- `remove`: テーブル型 (optional)
  - 系全体の並進・回転速度成分を取り除くことができます。`true`にすると取り除かれます。
  - `rescale`が`true`になっていた場合、全体の速度ベクトルをリスケールすることで速度を減算した分の運動エネルギーを補填します。
  - 省略した場合、全て`false`になります。
- `gammas`: テーブルの配列型
  - それぞれの粒子について、摩擦係数$$\gamma_i$$を指定します。
