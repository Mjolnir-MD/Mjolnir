+++
title = "BAOABLangevin"
weight = 1000
+++

# BAOABLangevin

ランジュバン方程式に従い、温度・体積・粒子数 (NVT) 一定のシミュレーションを行います。

{{<katex display>}}
m\frac{d^2 \bold{r}}{dt^2} = \bold{f}(\bold{r}) - m\gamma\bold{v} + \beta(t)
{{</katex>}}

以下の論文で提案された手法です。

- [Benedict Leimkuhler and Charles Matthews. Appl. Math. Res. Exp. (2013) 2013:1, pp. 34-56](https://doi.org/10.1093/amrx/abs010)
- [Benedict Leimkuhler and Charles Matthews. J. Chem. Phys. (2013) 138:17, 174102](https://doi.org/10.1063/1.4802990)

## 例

```toml
[simulator]
integrator.type = "BAOABLangevin"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## 入力

`delta_t`などの他のパラメータは[Simulator]({{<relref "/docs/reference/simulators">}})で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"BAOABLangevin"`です。
- `remove`: テーブル型 (optional)
  - 系全体の並進・回転速度成分を取り除くことができます。`true`にすると取り除かれます。
  - `rescale`が`true`になっていた場合、全体の速度ベクトルをリスケールすることで速度を減算した分の運動エネルギーを補填します。
  - 省略した場合、全て`false`になります。
- `gammas`: テーブルの配列型
  - 粒子の摩擦係数{{<katex>}}\gamma_i{{</katex>}}を指定します。

