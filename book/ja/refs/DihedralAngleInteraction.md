# DihedralAngleInteraction

`DihedralAngle`相互作用は、二面角にかかる相互作用です。

粒子`i`番目, `j`番目, `k`番目, `l`番目の粒子について、
`i`, `j`, `k`がなす面と`j`, `k`, `l`のなす面のなす角に対してかかります。

- [Harmonic](HarmonicPotential.md)
- [ClementiDihedral](ClementiDihedralPotential.md)
- [Gaussian](GaussianPotential.md)
- [FlexibleLocalDihedral](FlexibleLocalDihedral.md)

## 例

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Gaussian"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2, 3], ... }, # ポテンシャルによって要求されるパラメータは変化します。
    # 必要に応じて続きます...
]
```

## 入力

`DihedralAngle`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 結合角相互作用を使用する際は、`"BondAngle"`を指定します。
- `potential`: 文字列型
  - [`"Harmonic"`](HarmonicPotential.md): 調和振動子ポテンシャルを用います。
  - [`"ClementiDihedral"`](ClementiDihedralPotential.md): Clementi-Go力場に使用するポテンシャルです。
  - [`"Gaussian"`](GaussianPotential.md): ガウシアン型のポテンシャルを用います。
  - [`"FlexibleLocalDihedral"`](FlexibleLocalDihedral.md): AICG2+力場に使用するポテンシャルです。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - 他のパラメータは、ポテンシャルによって異なります。
