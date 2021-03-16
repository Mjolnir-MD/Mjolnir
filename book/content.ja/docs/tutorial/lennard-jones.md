+++
title  = "The Lennard-Jones fluid"
weight = 100
+++

# The Lennard-Jones fluid

簡単な例として、Lennard-Jonesポテンシャルで相互作用する粒子のシミュレーションをしてみましょう。
あまり時間がかかるシミュレーションにしたくないので、8x8x8の512粒子の系を考えることにします。
特に具体的な何かの原子をシミュレーションしたいというわけではないので、パラメータは適当に決めていきます。

このシミュレーションに必要なものは、大まかに

1. どのような方法でシミュレーションするか
2. 初期構造
3. ポテンシャルのパラメータ

です。`lennard-jones.toml`というファイルを作って、この順に実際に入力を書いていきましょう。

## `[files]` and `[units]`

具体的な設定を埋めていく前に、細々とした設定を先に済ませておきます。

まず、出力ファイルの名前とフォーマットを決めましょう。
名前はそのまま`lennard-jones`に、出力先はカレントディレクトリ、フォーマットは[xyz](https://en.wikipedia.org/wiki/XYZ_file_format)にします。

```toml
[files]
output.prefix = "lennard-jones"
output.path   = "./"
output.format = "xyz"
```

他に何が設定できるかは、[Referenceのfiles]({{<ref "/docs/reference/files">}})を見てください。

今後の入力で使う数値の単位系も決めないといけません。
ここでは、エネルギーの単位は`kcal/mol`、長さの単位は`angstrom`にします。

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

他に何が設定できるかは、[Referenceのunits]({{<ref "/docs/reference/units">}})を見てください。

## `[simulator]` and other general properties

さて、一般的な部分の設定は済ませたので、いよいよシミュレーションをどのように行うかを設定します。

まず、今から行うのは普通の分子動力学シミュレーションなので、方法は`MolecularDynamics`です。
浮動小数点数の型には倍精度を、境界条件には`Peirodic`を適用します。
それと、乱数シードを適当に決めておいてください。

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Periodic"
seed          = 123456789
```

あとは、時間発展をどうやって、どの程度行うかを決めます。

ここでは、ランジュバン方程式の時間積分法のうち、比較的最近提案されたBAOAB Langevin integratorを使うことにしましょう。
ここでは摩擦係数を粒子ごとに設定することができますが、今回は均質な系にするので、とりあえず全部1にしてしまいましょう。

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Periodic"
seed          = 123456789
delta_t       = 0.01
total_step    = 100_000
save_step     =     100
integrator.type = "BAOABLangevin"
integrator.gammas = [
    {index =   0, gamma = 1.00},
    {index =   1, gamma = 1.00},
    {index =   2, gamma = 1.00},
    # ...
    {index = 511, gamma = 1.00},
]
```

ちょっと大変ですが、エディタのマクロを使ったり、適当なスクリプトを書いてしまえばすぐです。

他に何が設定できるかは、[Referenceのsimulator]({{<ref "/docs/reference/simulators">}})を見てください。

## `[[systems]]`

続いて粒子の初期構造と、系全体のパラメータを設定しましょう。

まず、ランジュバン動力学で使う系の温度と、周期境界条件で使う境界の大きさを指定します。

```toml
[[systems]]
attributes.temperature = 300.0 # K
boundary_shape.lower  = [ 0.0,  0.0,  0.0]
boundary_shape.upper  = [16.0, 16.0, 16.0]
```

簡単に済ませたいので、粒子は格子状に並べることにします。
初期構造であまりにもエネルギーが高くなっては困るので、粒子の半径を1Åにして、距離2Åで並べていくことにしましょうか。

```toml
[[systems]]
attributes.temperature = 300.0 # K
boundary_shape.lower  = [ 0.0,  0.0,  0.0]
boundary_shape.upper  = [16.0, 16.0, 16.0]
particles = [
    {mass = 1.0, position = [ 1.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 3.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 5.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 7.000, 1.000, 1.000]},
    {mass = 1.0, position = [ 9.000, 1.000, 1.000]},
    {mass = 1.0, position = [11.000, 1.000, 1.000]},
    {mass = 1.0, position = [13.000, 1.000, 1.000]},
    {mass = 1.0, position = [15.000, 1.000, 1.000]},

    {mass = 1.0, position = [ 1.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 3.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 5.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 7.000, 3.000, 1.000]},
    {mass = 1.0, position = [ 9.000, 3.000, 1.000]},
    {mass = 1.0, position = [11.000, 3.000, 1.000]},
    {mass = 1.0, position = [13.000, 3.000, 1.000]},
    {mass = 1.0, position = [15.000, 3.000, 1.000]},
    # 続く...
]
```

他に何が設定できるかは、[ReferenceのSystems]({{<ref "/docs/reference/system">}})を見てください。

## `[[forcefields]]`

最後に、力場のパラメータを設定して終わりです。

Lennard-Jonesは非局所的な相互作用であり、かつ外力ではないので、`global`な相互作用に属します。
まずは`[[forcefields.global]]`を定義しましょう。
そしてそこに、二体間相互作用であることと、ポテンシャル関数としてLennard-Jonesを使うことを指定します。

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
```

続いて、計算の効率化のための空間分割法について設定します。
今回は、系が小さいのでCell Listでは速度が出ないかも知れません。
とりあえず、単純にリストを構築して一定時間使いまわすだけのVerletListにしておきます。

近距離ペアリストによる高速化では、カットオフ距離からマージンを取って少し離れた粒子までをリストに入れておき、粒子が移動した距離に応じてマージンを減らして、リストを再構築するまでの時間を自動的に計算します。
そのために、このマージンをどの程度にするかも決めなければなりません。
相互作用と系がどの程度密かに依存しますが、今回は適当に決めます。
これは計算結果に影響しない（計算が終わるまでの時間には影響します）ので、適当でも構いません。

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
spatial_partition.type = "VerletList"
spatial_partition.margin = 0.25
```

最後に、各粒子のパラメータを設定していきます。
二つの粒子の間にかかるポテンシャルのパラメータは、[Lorentz-Berthelot則](https://ja.wikipedia.org/wiki/%E3%83%AD%E3%83%BC%E3%83%AC%E3%83%B3%E3%83%84-%E3%83%99%E3%83%AB%E3%83%86%E3%83%AD%E5%89%87)によって決まります。

```toml
[[forcefields]]
[[forcefields.global]]
interaction = "Pair"
potential   = "LennardJones"
spatial_partition.type = "VerletList"
spatial_partition.margin = 0.25
parameters = [
    {index =   0, sigma = 2.0, epsilon = 1.0},
    {index =   1, sigma = 2.0, epsilon = 1.0},
    {index =   2, sigma = 2.0, epsilon = 1.0},
    # ...
    {index = 511, sigma = 2.0, epsilon = 1.0},
]

```

他に何が設定できるかは、[ReferenceのForceFields]({{<ref "/docs/reference/forcefields">}})を見てください。

## Simulation

お疲れ様でした！　これで入力ファイルは全部です。
このファイルをコンパイルしたmjolnirに渡せば、以下のようなログが出て、計算が始まります。
コンピュータによりますが、数分で終わるはずです。

```console
$ ./bin/mjolnir lennard-jones.toml
reading input file...
-- reading and parsing toml file `lennard-jones.toml` ...  successfully parsed.
-- the log file is `./lennard-jones.log`
-- mjolnir version v1.22.0
-- compiled using /usr/bin/g++-10
-- input_path is ./
-- expanding include files ...
-- done.
-- precision is double
-- Boundary Condition is Periodic. The shape is cuboid.
-- execute on single core
-- energy unit is [kcal/mol]
-- length unit is [angstrom]
-- Integrator is BAOABLangevin.
-- Simulator type is MolecularDynamics.
-- total step is 100000
-- save step is 100
-- checkpoint step is 100
-- reading 0th system ...
-- 512 particles are found.
-- output file prefix is `./lennard-jones`
-- output xyz format.
-- GlobalForceField (x1) found
-- reading 0th [[forcefields.global]]
-- Pair interaction found.
-- -- potential function is Lennard-Jones.
-- -- Spatial Partition is VerletList with relative margin = 0.25
-- No `ignore.group` is provided. All the groups are taken into account.
-- No `ignore.molecule` is provided. All the molecules are taken into account.
-- seed is 123456789
done.
initializing simulator...
-- generating velocity with T = 300...
-- done.
done.
start running simulation
 14.1%|███████                                           | 56.0 seconds remaining
```

計算が終われば、以下のようなファイルが出力されているはずです。

```console
$ ls lennard-jones*
lennard-jones.toml
lennard-jones.log
lennard-jones_velocity.xyz
lennard-jones_system.msg
lennard-jones_rng.msg
lennard-jones_position.xyz
lennard-jones.ene
```

`.msg`ファイルはシミュレーションを再開するためのリスタート用ファイルで、そのフォーマットは[MsgPack](https://msgpack.org/ja.html)です。

`.ene`ファイルはエネルギーなどの値が単純なテキストベースで書かれており、`gnuplot`などで簡単にプロットすることができます。

```console
$ head lennard-jones.ene
# unit of length : angstrom, unit of energy : kcal/mol
# timestep  GlobalPairLennardJones  kinetic_energy attribute:temperature
0                     -1704.786330     453.454993            300.000000
100                   -2579.996798     812.198804            300.000000
200                   -2787.935537     554.446648            300.000000
300                   -2909.224864     473.619251            300.000000
400                   -2913.189150     453.464065            300.000000
500                   -2964.234777     463.400579            300.000000
600                   -2988.270454     462.704726            300.000000
700                   -2960.111833     458.723132            300.000000
```

`lennard-jones_position.xyz`が位置のトラジェクトリです。
これをVMDなどに渡すことで、トラジェクトリを見ることができます。

## Conclusion

このチュートリアルは以上です。

入力ファイルの書き方を少し冗長だと思った人もいるかもしれません。
しかし、必要なパラメータを全て陽に指定できるということは、逆に簡単に好きなパラメータでシミュレーションを走らせることができるということでもあります。
このチュートリアルも、興味がある人はパラメータを一部変えてみて、例えば一つの粒子だけをとんでもなく大きくしたりして、何が起きるか観察して楽しんでみてください。
