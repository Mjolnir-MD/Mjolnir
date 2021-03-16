+++
title  = "A simple polymer model"
weight = 200
+++

# A simple polymer model

[前回]({{<relref "lennard-jones">}})に続いて、今度は局所的な相互作用の設定方法を、簡単なポリマーのシミュレーションを通して説明します。

今回は100個の粒子を調和振動子ポテンシャルで繋げてみます。
簡単のため、初めは繋がっていない粒子に非局所的な相互作用はないとします。

## `[files]` and `[units]`

前回と同様、細々とした設定を先に済ませてしまいましょう。

まずは、出力ファイルの名前とフォーマットを決めます。

```toml
[files]
output.prefix = "polymer-model"
output.path   = "./"
output.format = "xyz"
```

他に何が設定できるかは、[Referenceのfiles]({{<ref "/docs/reference/files">}})を見てください。

単位は、`kcal/mol`と`angstrom`にします。

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

他に何が設定できるかは、[Referenceのunits]({{<ref "/docs/reference/units">}})を見てください。

## `[simulator]`

ここもほぼ[前回]({{<relref "lennard-jones">}})と同じです。

ですが、今回は全ての粒子が繋がっていて、遠くに離れていく心配がないので、周期境界条件をなくしてみましょう。

```toml
[simulator]
type          = "MolecularDynamics"
precision     = "double"
boundary_type = "Unlimited" # No periodic boundary
seed          = 123456789
delta_t       = 0.01
total_step    = 1_000_000
save_step     =     1_000
integrator.type = "BAOABLangevin"
integrator.gammas = [
    {index =  0, gamma = 1.00},
    {index =  1, gamma = 1.00},
    {index =  2, gamma = 1.00},
    # ...
    {index = 99, gamma = 1.00},
]
```

他に何が設定できるかは、[ReferenceのSimulator]({{<ref "/docs/reference/simulators">}})を見てください。

## `[[systems]]`

[前回]({{<relref "lennard-jones">}})と違い周期境界条件を適用しないことにしたので、`boundary_shape`を指定する必要はありません。

今回は、最初は伸び切った直線のコンフォメーションからシミュレーションを始めます。

```toml
[[systems]]
attributes.temperature = 300.0 # K
particles = [
    {mass = 1.0, position = [ 0.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 1.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 2.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 3.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 4.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 5.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 6.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 7.000, 0.000, 0.000]},
    {mass = 1.0, position = [ 8.000, 0.000, 0.000]},
    # 続く...
    {mass = 1.0, position = [99.000, 0.000, 0.000]},
]
```

他に何が設定できるかは、[ReferenceのSystems]({{<ref "/docs/reference/system">}})を見てください。

## `[[forcefields]]`

最後に、力場のパラメータを設定します。

今回は結合長ポテンシャルを設定します。
これは局所的な（最初に指定したペアの間のみで働く）相互作用なので、定義するのは`[[forcefields.local]]`です。

相互作用の種類は結合長なので`BondLength`、ポテンシャル関数は調和振動子なので`Harmonic`です。

```toml
[[forcefields]]
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
```

`topology`というのは、あとで説明しますが、局所的な相互作用と非局所的な相互作用の間をうまく取り持つためのラベルです。
今回は必要ないので、適当な名前で埋めておきます。

さて、では隣り合う粒子の間にポテンシャルを適用していきましょう。

```toml
[[forcefields]]
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [ 0, 1], v0 = 1.0, k = 10.0},
    {indices = [ 1, 2], v0 = 1.0, k = 10.0},
    {indices = [ 2, 3], v0 = 1.0, k = 10.0},
    # ...
    {indices = [98,99], v0 = 1.0, k = 10.0},
]
```

これで力場のパラメータは全部です。

他に何が設定できるかは、[ReferenceのForceFields]({{<ref "/docs/reference/forcefields">}})を見てください。

## Simulation

お疲れ様でした！　これで入力ファイルは全部です。
このファイルをコンパイルしたmjolnirに渡せば、以下のようなログが出て、計算が始まります。

前回より粒子数が少ないので、すぐに終わるはずです。

```console
$ ./bin/mjolnir polymer-model.toml
reading input file...
-- reading and parsing toml file `polymer-model.toml` ...  successfully parsed.
-- the log file is `./polymer-model.log`
-- mjolnir version v1.22.0-dev (06d5ede6)
-- compiled using /usr/bin/g++-10
-- input_path is ./
-- expanding include files ...
-- done.
-- precision is double
-- Boundary Condition is Unlimited
-- execute on single core
-- energy unit is [kcal/mol]
-- length unit is [angstrom]
-- Integrator is BAOABLangevin.
-- Simulator type is MolecularDynamics.
-- total step is 1000000
-- save step is 1000
-- checkpoint step is 1000
-- reading 0th system ...
-- 100 particles are found.
-- output file prefix is `./polymer-model`
-- output xyz format.
-- LocalForceField (x1) found
-- reading 0th [[forcefields.local]]
-- Bond Length interaction found.
-- -- potential function is Harmonic.
-- -- 99 interactions are found.
-- seed is 123456789
done.
initializing simulator...
-- generating velocity with T = 300...
-- done.
done.
start running simulation
  8.2%|████                                              | 10.5 seconds remaining
```

計算が終われば、以下のようなファイルが出力されているはずです。

```console
$ ls polymer-model*
polymer-model.toml
polymer-model.ene
polymer-model.log
polymer-model_rng.msg
polymer-model_system.msg
polymer-model_position.xyz
polymer-model_velocity.xyz
```

`.msg`ファイルはシミュレーションを再開するためのリスタート用ファイルで、そのフォーマットは[MsgPack](https://msgpack.org/ja.html)です。

`.ene`ファイルはエネルギーなどの値が単純なテキストベースで書かれており、`gnuplot`などで簡単にプロットすることができます。

```console
$ head polymer-model.ene
# unit of length : angstrom, unit of energy : kcal/mol
# timestep  BondLength:Harmonic  kinetic_energy attribute:temperature
0                      0.000000      96.956625            300.000000
1000                  56.910239      86.124451            300.000000
2000                  51.497857      85.960034            300.000000
3000                  44.394914      90.433666            300.000000
4000                  50.774066      89.091169            300.000000
5000                  46.749105     104.594276            300.000000
6000                  34.130104      86.045558            300.000000
7000                  39.505721      84.893861            300.000000
```

`polymer-model_position.xyz`が位置のトラジェクトリです。
これをVMDなどに渡すことで、トラジェクトリを見ることができます。

## Conclusion

このチュートリアルは以上です。

`BondLength`を設定する二つの粒子の番号を直接指定しました。
これはなにも連続した二粒子である必要はなく、好きな粒子のペアに相互作用をかけることができます。
最初と最後の粒子を繋いで環状にすることも、一つ飛ばしで繋いで角度をつけることもできます。

局所的な相互作用は他にも、結合角、二面角などが知られています。
Mjolnirはこれらを同様の形でサポートしているので、興味のある組み合わせを試してみてください。

それらの設定方法は、[ReferenceのLocalForceField]({{<ref "docs/reference/forcefields/local">}})を見てください。
