# 3SPN2BaseBaseInteraction

`3SPN2BaseBase`相互作用は、3SPN2力場で使用される`BasePairing`と`CrossStacking`の
両方を同時に計算する相互作用です。

利用可能なポテンシャルは以下のものがあります。

- `"3SPN2"`
  - Hinckley et al., (2013) JCP
- `"3SPN2C"`
  - Freeman et al., (2014) JCP

## 例

```toml
[[forcefields.global]]
interaction = "3SPN2BaseBase"
potential   = "3SPN2"
ignore.particles_within.nucleotide = 3
spatial_partition = {type = "CellList", margin = 0.2}
parameters  = [
# `nucleotide` index starts from 5' and ends at 3'.
{strand = 0, nucleotide =  0,          S =   0, B =   1, offset = 100, Base = "A"},
{strand = 0, nucleotide =  1, P =   2, S =   3, B =   4, offset = 100, Base = "T"},
{strand = 0, nucleotide =  2, P =   5, S =   6, B =   7, offset = 100, Base = "C"},
# ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の種類を設定します。`"3SPN2BaseBase"`です。
- `ignore`: テーブル型
  - 相互作用を無視する条件を記述します。
  - 3SPN2力場では3つのヌクレオチドが間にないと塩基対を形成しないので、`BaseStacking`の
    トポロジーとして設定した結合(ここでは`"nucleotide"`)を3つ分無視して下さい。
- `spatial_partition`: テーブル型
  - 計算速度の向上のための空間分割の方法を指定します。
- `potential`: 文字列型
  - ポテンシャルの種類を設定します。
  - `3SPN2` : 3SPN.2ポテンシャルのパラメータを使います。
  - `3SPN2C`: 3SPN.2Cポテンシャルのパラメータを使います。
- `parameters`:テーブルの配列型
  - 力場が適用される粒子ごとにパラメータを与えます。ここでは、
    - ヌクレオチドが所属するストランドのインデックス
    - ヌクレオチドの通し番号
    - リン酸、糖、塩基に対応する粒子の番号、
    - 塩基の種類
  - を指定します。
  - `index`: 整数型
    - どの粒子に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - 各粒子の粒子番号に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。

### `ignore`

相互作用を無視する条件を指定します。
[`GlobalPair`](GlobalPairInteraction.md)でのそれと同じものが使えます。

### `spatial_partition`

高速化のための空間分割の方法を定義します。
[`GlobalPair`](GlobalPairInteraction.md)でのそれと同じものが使えます。
