# SimulatedAnnealing

焼きなまし（Simulated Annealing）法のためのシミュレータです。

一つの[`System`](System.md)、一つの[`ForceField`](ForceField.md)を指定して、
それを使って焼きなまし法を行います。

## 例

```toml
[simulator]
type           = "SimulatedAnnealing"
boundary_type  = "Unlimited"
precision      = "double"
parallelism    = "OpenMP" # optional
delta_t        = 0.1
total_step     = 50_000
save_step      = 100
each_step      = 100
schedule.type  = "linear"
schedule.begin = 300
schedule.end   = 150
[simulator.integrator]
type = "UnderdampedLangevin"
# ...
```

## 入力

以下のパラメータを取ります。

- `type`: 文字列型
  - シミュレータの種類を指定します。このシミュレータを使う場合、`"SimulatedAnnealing"`です。
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
- `each_step`: 整数型
  - 何ステップおきに温度を変更するか指定します。
- `schedule`: テーブル型
  - どのように温度を変更するか指定します。
- `integrator`: テーブル型
  - 時間積分の方法を指定します。積分方法によって必要なパラメータが異なります。

### `schedule`

`simulator.schedule`テーブルは、以下のパラメータを取ります。

- `type`: 文字列型
  - 温度の変化させ方を指定します。以下の種類がサポートされています。
  - `"linear"`: 線形に変化させます。
- `begin`: 浮動小数点数型
  - `"linear"`な温度変化において、時刻0での温度を指定します。単位は[K]です。
- `end`: 浮動小数点数型
  - `"linear"`な温度変化において、最終時刻での温度を指定します。単位は[K]です。

### `integrator`

- `type`: 文字列型
  - 積分方法の種類を指定します。`MolecularDynamics`とは異なり、温度一定の手法しか使用できません。
    - ["UnderdampedLangevin"](UnderdampedLangevinIntegrator.md)
      - NVT一定のランジュバン方程式に基づくシミュレーションを行います。
