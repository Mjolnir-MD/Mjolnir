+++
title = "VelocityVerlet"
weight = 5000
+++

# VelocityVerlet

ニュートンの運動方程式に従い、エネルギー・体積・粒子数一定のシミュレーションを行います。

## 例

```toml
[simulator]
integrator.type = "VelocityVerlet"
integrator.remove.translation = true
integrator.remove.rotation    = true
integrator.remove.rescale     = true
```

## 入力

`delta_t`などの他のパラメータは[`[simulator]`]({{<relref "/docs/reference/simulators">}})で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"VelocityVerlet"`です。
- `remove`: テーブル型 (optional)
  - `translation`: 論理値型
    - `true`の場合、毎ステップ、系全体の並進速度成分を取り除きます。
  - `rotation`: 論理値型
    - `true`の場合、毎ステップ、系全体の回転速度成分を取り除きます。
  - `rescale`: 論理値型
    - `true`になっていた場合、全体の速度ベクトルをリスケールすることで速度を減算した分の運動エネルギーを補填します。
  - 省略した場合、全て`false`になります。
