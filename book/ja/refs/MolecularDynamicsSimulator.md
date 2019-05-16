# MolecularDynamics

通常の分子動力学シミュレーションを行うシミュレータです。

一つの[`System`](System.md)、一つの[`ForceField`](ForceField.md)を指定して、
それを使って分子動力学シミュレーションを行います。

## 例

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "PeriodicCuboid"
precision     = "double"
parallelism   = "OpenMP" # optional
delta_t       = 0.1
total_step    = 50_000
save_step     = 100
integrator.type = "VelocityVerlet"
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
- `delta_t`: 浮動小数点数型
  - シミュレーションの時間刻みを指定します。
  - 時間の単位は[`units`](Units.md)で指定した単位系に依存します。
- `total_step`: 整数型
  - 実行するステップ数を指定します。
- `save_step`: 整数型
  - 何ステップおきに状態を出力するか指定します。
- `integrator`: テーブル型
  - 時間積分の方法を指定します。積分方法によって必要なパラメータが異なります。

### `integrator`

- `type`: 文字列型
  - 積分方法の種類を指定します。現在、以下の積分方法が利用可能です。
    - ["UnderdampedLangevin"](UnderdampedLangevinIntegrator.md)
      - NVT一定のランジュバン方程式に基づくシミュレーションを行います。
    - ["VelocityVerlet"](VelocityVerletIntegrator.md)
      - エネルギー一定のニュートン力学に基づくシミュレーションを行います。
