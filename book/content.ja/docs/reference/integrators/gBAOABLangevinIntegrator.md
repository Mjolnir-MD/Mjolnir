+++
title = "g-BAOABLangevin"
weight = 2000
+++

# g-BAOABLangevin

ランジュバン方程式に従い、温度・体積・粒子数 (NVT) 一定のシミュレーションを行います。

`BAOABLangevin`と異なり、結合長に対する拘束条件を取り扱うことができます。

以下の論文で提案された手法です。

- [Leimkuhler B, Matthews C. Proc. R. Soc. A. (2016)](https://doi.org/10.1098/rspa.2016.0138)

## 例

```toml
[simulator]
# ...
integrator.type = "g-BAOABLangevin"
integrator.gammas = [
    {index = 0, gamma = 1.0},
    {index = 1, gamma = 1.0},
    # ...
]
```

## 入力

`delta_t`などの他のパラメータは[Simulator]({{<relref "/docs/reference/simulators">}})で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"g-BAOABLangevin"`です。
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

## Remarks

This feature is developed by contributor, [@yutakasi634](https://github.com/yutakasi634).
