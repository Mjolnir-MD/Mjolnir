+++
title = "ForceFields"
weight = 600
bookCollapseSection = true
+++

# `[[forcefields]]`

用いる力場を設定します。
パラメータの単位は[`[units]`]({{<relref "/docs/reference/units">}})での設定に準じます。

複数の力場を扱うシミュレーションのために、テーブルの配列として定義されています。
通常のシミュレーションでは、一つしか定義されません。

## env

`[[forcefields]]`内では、特殊な値`env`を使用することができます。

`env`はテーブルで、ここに値を定義することで、同じ相互作用テーブルの中で値を参照することができます。

```toml
[[forcefields.local]]
env.pi = 3.1416 # これによって3.1416の代わりに文字列"pi"を使用できるようになります。
parameters = [
    {indices = [0, 1], k = 100.0, v0 = "pi"},
    {indices = [1, 2], k = 100.0, v0 = "pi"},
    {indices = [2, 3], k = 100.0, v0 = "pi"},
]
```

この機能は[ファイルのインクルード]({{<relref "/docs/reference/_index.md#%E3%83%95%E3%82%A1%E3%82%A4%E3%83%AB%E3%81%AE%E3%82%A4%E3%83%B3%E3%82%AF%E3%83%AB%E3%83%BC%E3%83%89">}})と組み合わせて使用することができます。

## [LocalForceFiled]({{<relref "local">}})

決まった粒子の間のみにかかる相互作用です。結合長、結合角、二面角などが該当します。

## [GlobalForceFiled]({{<relref "global">}})

対応する粒子の全てのペアにかかる相互作用です。ペア相互作用、例えば静電相互作用などが該当します。

## [ExternalForceFiled]({{<relref "external">}})

外力による相互作用です。空間中の点に粒子を束縛する相互作用や、壁を設定する相互作用が該当します。

## [ConstraintForceField]({{<relref "constraint">}})

拘束条件です。特定の粒子間の距離を一定に保ちます。

## [MultipleBasinForceField]({{<relref "MultipleBasinForceField.md">}})

特殊なメタ力場で、複数の力場をスムースに結合することができます。
