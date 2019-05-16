# Simulator

`[simulator]`テーブルでは、シミュレーション全体の設定と時間積分のアルゴリズムを設定します。

全シミュレータに共通して、シミュレーションに用いる浮動小数点数の精度と、境界条件の種類を指定する必要があります。

## [MolecularDynamicsSimulator](MolecularDynamicsSimulator.md)

ごく通常のMDシミュレーションを行うためのシミュレータです。

## [SimulatedAnnealingSimulator](SimulatedAnnealingSimulator.md)

温度を変えながらシミュレーションすることによって効率よく最安定点を探す、Simulated Annealing法のシミュレータです。

## [SteepestDescentSimulator](SteepestDescentSimulator.md)

力場の勾配によって最安定点を探す、Steepest Descent法のシミュレータです。

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
