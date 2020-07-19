+++
title = "SwitchingForceField"
weight = 4000
+++

# SwitchingForceField

あらかじめ決めた時間ステップに使用する力場を変更するシミュレータです。

[system]({{<relref "/docs/reference/system">}})を一つ、[forcefield]({{<relref "/docs/reference/forcefields">}})を複数個（必要なだけ）要求します。

## Example

```toml
[simulator]
type          = "SwitchingForceField"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 3_000_000
save_step     = 100
seed          = 2374

integrator.type = "BAOABLangevin"
integrator.parameters = [
{index =  0, gamma = 1.00},
{index =  1, gamma = 1.00},
]

schedule      = [
    {until = 1_000_000, forcefield = "close"},
    {until = 2_000_000, forcefield = "open"},
    {until = 3_000_000, forcefield = "close"},
]

[[forcefields]]
name = "close"
[[forcefields.local]]
# ...

[[forcefields]]
name = "open"
[[forcefields.local]]
# ...
```

## Input Reference


- `type`: 文字列型
  - シミュレータの種類を指定します。このシミュレータを使う場合、`"SwitchingForceField"`です。
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
- `seed`: 整数型
  - 乱数生成器の初期化に用いるシードを設定します。
- `delta_t`: 浮動小数点数型
  - シミュレーションの時間刻みを指定します。
  - 時間の単位は[`units`](Units.md)で指定した単位系に依存します。
- `total_step`: 整数型
  - 実行するステップ数を指定します。
- `save_step`: 整数型
  - 何ステップおきに状態を出力するか指定します。
- `integrator`: テーブル型
  - 時間積分の方法を指定します。積分方法によって必要なパラメータが異なります。
  - 以下のアルゴリズムが使用可能です。
  - ["BAOABLangevin"]({{<relref "/docs/reference/integrators/BAOABLangevinIntegrator.md">}})
  - ["g-BAOABLangevin"]({{<relref "/docs/reference/integrators/gBAOABLangevinIntegrator.md">}})
  - ["UnderdampedLangevin"]({{<relref "/docs/reference/integrators/UnderdampedLangevinIntegrator.md">}})
  - ["VelocityVerlet"]({{<relref "/docs/reference/integrators/VelocityVerletIntegrator.md">}})
  - 参考：[`integrators`]({{<relref "/docs/reference/integrators">}}).
- `schedule`: テーブルの配列型
  - 使用する力場のスケジュールを決めます。以下で詳細を説明します。

### `schedule`

テーブルの配列で、各テーブルが以下のパラメータを持ちます。

- `until`: 整数型
  - この時間ステップになるまで、指定された力場を使用します。
  - この時間ステップ以降は、次の力場を使用します。
- `forcefield`: 文字列型
  - 使用する力場の名前を指定します。

`forcefield`には、`[[forcefields]]`直下の`name`で指定した名前を使って下さい。

```toml
[[forcefields]]
name = "close"

[[forcefields.local]]
interaction = "BondLength"
potential = "Harmonic"
# ここは`close`の一員になります

[[forcefields.local]]
interaction = "BondAngle"
potential = "Harmonic"
# ここも`close`の一員になります
# ...

[[forcefields]]
name = "open"

[[forcefields.local]]
# これ以降は`open`の一員になります
# ...
```
