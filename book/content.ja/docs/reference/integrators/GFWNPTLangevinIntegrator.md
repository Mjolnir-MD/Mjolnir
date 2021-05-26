+++
title = "GFWNPTLangevin"
weight = 4000
+++

# GFWNPTLangevin

ランジュバン方程式に従い、温度・圧力・粒子数 (NPT) 一定のシミュレーションを行います。

{{<katex display>}}
m\frac{d^2 \bold{r}}{dt^2} = \bold{f}(\bold{r}) - m\gamma\bold{v} + \beta(t)
{{</katex>}}

以下の論文で提案された手法です。

- [Xingyu Gao, Jun Fang, and Han Wang. J. Chem. Phys. (2016) 144, 124113](https://doi.org/10.1063/1.4944909)

## 例

```toml
[simulator]
# ...
integrator.type = "GFWNPTLangevin"
integrator.chi  = 0.0
integrator.cell_mass  = [1e3, 1e3, 1e3]
integrator.cell_gamma = [0.1, 0.1, 0.1]
integrator.gammas = [
    {index = 0, gamma = 0.1},
    {index = 1, gamma = 0.1},
    # ...
]
```

## 入力

`delta_t`などの他のパラメータは[Simulator]({{<relref "/docs/reference/simulators">}})で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"BAOABLangevin"`です。
- `chi`: 浮動小数点数
  - パラメータ{{<katex>}} chi {{</katex>}}の値です。
- `cell_mass`: 浮動小数点数の配列型
  - 箱の運動に対応する質量です。
- `cell_gamma`: 浮動小数点数の配列型
  - 箱の運動に対応する摩擦係数です。
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
