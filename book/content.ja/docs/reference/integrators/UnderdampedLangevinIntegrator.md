+++
title = "UnderdampedLangevin"
weight = 4000
+++

# UnderdampedLangevin

ランジュバン方程式に従い、温度・体積・粒子数一定のシミュレーションを行います。

{{<katex display>}}
m\frac{d^2 \bold{r}}{dt^2} = \bold{f}(\bold{r}) - m\gamma\bold{v} + \beta(t)
{{</katex>}}


以下の論文で提案された手法です。

- J. D. Honeycutt and D. Thirumalai, (1992) Biopolymers
- Z. Guo and D. Thirumalai, (1995) Biopolymers.

## 例

```toml
[simulator]
integrator.type = "UnderdampedLangevin"
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
- `gammas`: テーブルの配列型
  - 粒子の摩擦係数{{<katex>}}\gamma_i{{</katex>}}を指定します。
- `remove`: テーブル型 (optional)
  - `translation`: 論理値型
    - `true`の場合、毎ステップ、系全体の並進速度成分を取り除きます。
  - `rotation`: 論理値型
    - `true`の場合、毎ステップ、系全体の回転速度成分を取り除きます。
  - `rescale`: 論理値型
    - `true`になっていた場合、全体の速度ベクトルをリスケールすることで速度を減算した分の運動エネルギーを補填します。
  - 省略した場合、全て`false`になります。
