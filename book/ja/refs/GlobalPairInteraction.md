# GlobalPair

`GlobalPair`相互作用は、パラメータを設定された粒子群の間の全てのペアについて働く相互作用です。

利用可能なポテンシャルには以下のようなものがあります。

- [ExcludedVolume](ExcludedVolumePotential.md)
- [LennardJones](LennardJonesPotential.md)
- [UniformLennardJones](UniformLennardJonesPotential.md)
- [DebyeHuckel](DebyeHuckelPotential.md)

## 例

```toml
[[forcefields.global]]
interaction = "Pair"
potential   = "ExcludedVolume"
ignore.molecule = "Nothing"
ignore.particles_within.bond    = 3
ignore.particles_within.contact = 1
spatial_partition = {type = "CellList", margin = 0.2}
parameters = [
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の種類を設定します。ペア相互作用を使用する場合は、`"Pair"`です。
- `potential`: 文字列型
  - ポテンシャルの種類を設定します。
  - [`"ExcludedVolume"`](ExcludedVolumePotential.md): 排除体積ポテンシャルを用います。
  - [`"LennardJones"`](LennardJonesPotential.md): レナード・ジョーンズポテンシャルを用います。
  - [`"UniformLennardJones"`](UniformLennardJonesPotential.md): パラメータが粒子によらず一定なレナード・ジョーンズポテンシャルを用います。
  - [`"DebyeHuckel"`](DebyeHuckelPotential.md): 溶媒による遮蔽効果を考慮した静電相互作用を用います。
- `ignore`: テーブル型
  - 相互作用を無視する条件を記述します。
- `spatial_partition`: テーブル型
  - 計算速度の向上のための空間分割の方法を指定します。

### `ignore`

相互作用を無視する条件を指定します。

特定の`LocalInteraction`が定義されている粒子のペアに対して無視する方法と、
同一分子内、または分子間の相互作用を無視する方法があります。

これらのトポロジーをどのように定義するかの詳細については、[`Topology`](Topology.md)を参照して下さい。

- `ignore.molecule`: 文字列型
  - 分子内、分子間の相互作用を無視するかどうかを設定します。
  - `"Nothing"`: 全ての粒子間の相互作用を含めます。
  - `"Self"`: 同一分子内の相互作用を無視します。別の分子間の相互作用のみが考慮されます。
  - `"Others"`: 別の分子間の相互作用を無視します。同一分子内の相互作用のみが考慮されます。
- `ignore.particles_within`: テーブル型
  - `LocalInteraction`が定義されているペアの間の相互作用を無視します。
  - トポロジーの名称と、その相互作用をいくつまで辿って無視するかの本数を指定します。
  - 例えば、`bond = 3`と指定した場合、トポロジーが`bond`と定義された相互作用が3つまでで繋がっているペアは無視されます。

### `spatial_partition`

高速化のための空間分割の方法を定義します。

- `type`: 文字列型
  - 空間分割の方法を指定します。以下の種類があります。
  - `"CellList"`: セルリストを用いてリストを作成します。
  - `"VerletList"`: 全粒子間の距離を計算してリストを作成します。
- `margin`: 浮動小数点数型
  - リストを作成する際のマージンの長さを、カットオフとの相対長さで指定します。
  - 粒子がこのマージン距離だけ動く毎にリストは更新されます。最適な長さは状況によって異なります。
  - 長くとると更新頻度が減りますが、リストが大きくなります。
  - 短くなるとリストが小さくなりますが、更新頻度が上がります。

