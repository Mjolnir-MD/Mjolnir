# BondAngleInteraction

`BondAngle`相互作用は、結合角にかかる相互作用です。

特定の粒子、`i`番目、`j`番目、`k`番目の粒子がなす角`ijk`にかかる角度です。
利用可能なポテンシャルには以下のようなものがあります。

- [Harmonic](HarmonicPotential.md)
- [Gaussian](GaussianPotential.md)
- [FlexibleLocalAngle](FlexibleLocalAngle.md)

## 例

```toml
[[forcefields.local]]
interaction = "BondAngle"
potential   = "Harmonic"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2], offset = 100, ... }, # ポテンシャルによって要求されるパラメータは変化します。
    # 必要に応じて続きます...
]
```

## 入力

`BondAngle`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 結合角相互作用を使用する際は、`"BondAngle"`を指定します。
- `potential`: 文字列型
  - [`"Harmonic"`](HarmonicPotential.md): 調和振動子ポテンシャルを用います。
  - [`"Gaussian"`](GaussianPotential.md): ガウシアン型のポテンシャルを用います。
  - [`"FlexibleLocalAngle"`](FlexibleLocalAngle.md): AICG2+力場に使用するポテンシャルです。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。
