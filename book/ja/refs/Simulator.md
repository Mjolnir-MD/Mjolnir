# Simulator

`[simulator]`テーブルでは、シミュレーション全体の設定と時間積分のアルゴリズムを設定します。

全シミュレータに共通して、シミュレーションに用いる浮動小数点数の精度と、境界条件の種類を指定する必要があります。

また、並列化を行う際もこの時点で設定を行います。

## Common part

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "PeriodicCuboid"
precision     = "double"
parallelism   = "OpenMP" # optional
```

全てのシミュレータは、共通して上記のパラメータを取ります。

- `type`: 文字列型
  - シミュレータの種類を指定します。個別のページを参照して下さい。
- `boundary_type`: 文字列型
  - 境界条件の種類を選択します。具体的な形状は[`System`](System.md)で指定します。
  - `"Unlimited"`: 境界条件を設定しません。シミュレーションボックスは無限大の大きさになります。
  - `"PeriodicCuboid"`: 直方体型の周期境界条件を指定します。
- `precition`: 文字列型
  - シミュレーションに用いる浮動小数点数型の種類を指定します。
  - `"float"`: 32bit浮動小数点数を使用します。
  - `"double"`: 64bit浮動小数点数を使用します。
- `parallelism`: 文字列型(省略可)
  - 並列化する際の実装を選択します。省略した際は、シングルコアで実行されます。
  - `"OpenMP"`: OpenMPを使った実装を使用します。
  - `"sequencial"`: 並列化を行いません。省略した場合はこれが選択されます。

通常、1コアで計算する場合は並列化オプションを使用しない方が高速です。

## [MolecularDynamicsSimulator](MolecularDynamicsSimulator.md)

ごく通常のMDシミュレーションを行うためのシミュレータです。

## [SimulatedAnnealingSimulator](SimulatedAnnealingSimulator.md)

温度を変えながらシミュレーションすることによって効率よく最安定点を探す、Simulated Annealing法のシミュレータです。

## [SteepestDescentSimulator](SteepestDescentSimulator.md)

力場の勾配によって最安定点を探す、Steepest Descent法のシミュレータです。

{% hint style='working' %}
Version 1.1.0 現在、`SteepestDescentSimulator`は並列化に対応していません。
{% endhint %}

## [SwitchingForceFieldSimulator](SteepestDescentSimulator.md)

あらかじめ決めた時間ステップに使用する力場を変更するシミュレータです。

## ファイル分割の方法

`[simulator]`テーブルは、メインの入力ファイルと事なるファイルに分割することが可能です。
以下のようにします。

```toml
# メインの入力ファイル
[files]
input.path = "./input/"

[simulator]
file_name = "simulator.toml"
```

`files.input.path`と`file_name`を結合した名前のファイルが読み込まれます。

```toml
# ./input/simulator.toml
[simulator]
type = "MolecularDynamics"
# ...
```
