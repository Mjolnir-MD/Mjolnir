# ForceField

用いる力場を設定します。
パラメータの単位は`[units]`での設定に準じます。

複数の力場を扱うシミュレーションのために、テーブルの配列として定義されています。
通常のシミュレーションでは、一つしか定義されません。

メインの入力ファイルとは別のファイルとして設定可能です。

## [LocalForceFiled](LocalForceField.md)

決まった粒子の間のみにかかる相互作用です。結合長、結合角、二面角などが該当します。

## [GlobalForceFiled](GlobalForceField.md)

対応する粒子の全てのペアにかかる相互作用です。ペア相互作用、例えば静電相互作用などが該当します。

## [ExternalForceFiled](ExternalForceFiled.md)

外力による相互作用です。空間中の点に粒子を束縛する相互作用や、壁を設定する相互作用が該当します。

## ファイル分割の方法

`[[forcefields]]`テーブルは、メインの入力ファイルと事なるファイルに分割することが可能です。
以下のようにします。

```toml
# メインの入力ファイル
[files]
input.path = "./input/"

[[forcefields]]
file_name = "forcefield.toml"
```

```toml
# ./input/forcefield.toml
[[forcefields]]
# ...
```

`files.input.path`と`file_name`を結合した名前のファイルが読み込まれます。

----

あるいは、もう少し細かく、相互作用ごとに分割することも可能です。

```toml
[[forcefields.local]]
file_name = "bond_length.toml"
[[forcefields.global]]
file_name = "electrostatic.toml"
```

```toml
# ./input/bond_length.toml
[[local]]
# ...
```

```toml
# ./input/electrostatic.toml
[[global]]]
# ...
```

ここで、テーブルの階層を合わせるため、`forcefields.local`ではなく単に`local`が
キーとして使用されています。
