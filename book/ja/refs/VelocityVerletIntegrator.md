# VelocityVerlet

ニュートンの運動方程式に従い、エネルギー・体積・粒子数一定のシミュレーションを行います。

## 例

```toml
[simulator]
integrator.type = "VelocityVerlet"
```

## 入力

`delta_t`などの他のパラメータは[Simulator](Simulator.md)で設定します。

- `type`: 文字列型
  - [Integrator](Integrator.md)の種類を指定します。`"VelocityVerlet"`です。
