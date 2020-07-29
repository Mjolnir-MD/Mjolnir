+++
title  = "Simulators"
weight = 300
bookCollapseSection = true
+++

# `[simulator]`

`[simulator]`テーブルでは、シミュレーション全体の設定と時間積分のアルゴリズムを設定します。

全シミュレータに共通して、シミュレーションに用いる浮動小数点数の精度と、境界条件の種類を指定する必要があります。

また、並列化を行う際もこの時点で設定を行います。

## Input Reference

### Common Part

以下のフィールドは全てのシミュレータに共通です。

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Periodic"
precision     = "double"
parallelism   = "OpenMP" # optional
```

- `type`: 文字列型
  - シミュレータの種類を指定します。例えば、普通のシミュレーションの他に、焼きなまし法やエネルギー最小化などが選べます。
  - 選択できるシミュレータは"Available Simulators"を参照してください。
- `boundary_type`: 文字列型
  - 境界条件の種類を選択します。具体的な形状は[`System`]({{<relref "/docs/reference/system">}})で指定します。
  - `"Unlimited"`: 境界条件を設定しません。シミュレーションボックスは無限大の大きさになります。
  - `"Periodic"`: 直方体型の周期境界条件を指定します。
  - `"PeriodicCuboid"`: 直方体型の周期境界条件を指定します。上と同じです。
- `precition`: 文字列型
  - シミュレーションに用いる浮動小数点数型の種類を指定します。
  - `"float"`: 32bit浮動小数点数を使用します。
  - `"double"`: 64bit浮動小数点数を使用します。
- `parallelism`: 文字列型(省略可)
  - 並列化する際の実装を選択します。省略した際は、シングルコアで実行されます。
  - `"OpenMP"`: OpenMPを使った実装を使用します。
  - `"sequencial"`: 並列化を行いません。省略した場合はこれが選択されます。
- `forcefield`: テーブル型 (省略可)
  - 特殊な場合のためのフィールドです。
  - 使用例は、[MultipleBasinForceField]({{<relref "/docs/reference/forcefields/MultipleBasinForceField.md">}})を参照してください。

### Available Simulators

- [MolecularDynamics]({{<relref "MolecularDynamicsSimulator.md">}})
  - ごく通常の分子動力学シミュレーションを行います。
- [SimulatedAnnealing]({{<relref "SimulatedAnnealingSimulator.md">}})
  - 焼きなまし法を行います。
- [SteepestDescent]({{<relref "SteepestDescentSimulator.md">}})
  - 最急降下法によってエネルギー極小の構造を探します。
- [SwitchingForceField]({{<relref "SwitchingForceFieldSimulator.md">}})
  - 指定した時刻に力場を変更するシミュレーションを行います。
- [EnergyCalculation]({{<relref "EnergyCalculationSimulator.md">}})
  - 与えられたトラジェクトリの各フレームでのエネルギーを、指定された力場で計算します。
