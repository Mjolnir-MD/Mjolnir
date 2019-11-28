# EnergyCalculation

トラジェクトリファイルからエネルギーを計算するシミュレータです。

[`ForceField`](ForceField.md)を指定して、トラジェクトリファイルのそれぞれの
スナップショットでのエネルギーを計算します。

このシミュレータに[`Integrator`](Integrator.md)は必要ありません。
また、[`System`](System.md)で指定された座標はトラジェクトリファイルの座標で上書きされます。

## 例

```toml
[simulator]
type          = "EnergyCalculation"
boundary_type = "PeriodicCuboid"
precision     = "double"
parallelism   = "OpenMP" # optional
file          = "example_position.dcd"
```

## 入力

以下のパラメータを取ります。

- `type`: 文字列型
  - シミュレータの種類を指定します。このシミュレータを使う場合、`"MolecularDynamics"`です。
- `boundary_type`: 文字列型
  - 境界条件の種類を指定します。具体的な大きさは[`System`](System.md)で指定します。
  - `"Unlimited"`: 境界条件を設定しません。シミュレーションボックスは無限大の大きさになります。
  - `"PeriodicCuboid"`: 直方体型の周期境界条件を指定します。
- `precision`: 文字列型
  - シミュレーションに用いる浮動小数点数型の種類を指定します。
  - `"float"`: 32bit浮動小数点数を使用します。
  - `"double"`: 64bit浮動小数点数を使用します。
- `parallelism`: 文字列型(省略可)
  - 並列化する際の実装を選択します。
  - `"OpenMP"`: OpenMPを使った実装を使用します。
  - `"sequencial"`: 並列化を行いません。省略した場合はこれが選択されます。
- `file`: 文字列型
  - エネルギーを計算するトラジェクトリファイルを指定します。
