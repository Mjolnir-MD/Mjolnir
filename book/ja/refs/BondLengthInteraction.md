# BondLengthInteraction

`BondLength`相互作用は、その名の通り結合長にかかる相互作用です。

特定の粒子、`i`番目と`j`番目の粒子の間の距離に対してかかる相互作用で、
利用可能なポテンシャルには以下のようなものがあります。

- [Harmonic](HarmonicPotential.md)
- [Gaussian](GaussianPotential.md)
- [GoContact](GoContactPotential.md)
- [WormLikeChain](WormLikeChainPotential.md)

## 例

```toml
[[forcefields.local]]
interaction = "BondLength"
potential   = "Harmonic"
topology    = "bond"
parameters  = [
    {indices = [0, 1], ... }, # ポテンシャルによって要求されるパラメータは変化します。
    # 必要に応じて続きます...
]
```

## 入力

`BondLength`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 結合長相互作用を使用する際は、`"BondLength"`を指定します。
- `potential`: 文字列型
  - [`"Harmonic"`](HarmonicPotential.md): 調和振動子ポテンシャルを用います。
  - [`"Gaussian"`](GaussianPotential.md): ガウシアン型のポテンシャルを用います。
  - [`"GoContact"`](GoContactPotential.md): L-J(10,12)型のポテンシャルを用います。
  - [`"WormLikeChain"`](WormLikeChainPotential.md): Worm-like chainモデルを元にしたポテンシャルを用います。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - 他のパラメータは、ポテンシャルによって異なります。
