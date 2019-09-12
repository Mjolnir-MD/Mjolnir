# ContactInteraction

`Contact`相互作用は、粒子間距離にかかる相互作用です。
`BondLength`相互作用と基本的に同じですが、シミュレーション中に力がかからない距離まで
粒子が移動する可能性を考慮し、力がかからなくなったら力の計算をスキップするようになっています。

特定の粒子、`i`番目と`j`番目の粒子の間の距離に対してかかる相互作用で、
利用可能なポテンシャルには以下のようなものがあります。

- [Gaussian](GaussianPotential.md)
- [GoContact](GoContactPotential.md)

## 例

```toml
[[forcefields.local]]
interaction = "Contact"
potential   = "GoContact"
topology    = "bond"
margin      = 0.5 # relative length to longest cutoff
parameters  = [
    {indices = [0, 1], ... }, # ポテンシャルによって要求されるパラメータは変化します。
    # 必要に応じて続きます...
]
```

## 入力

`Contact`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 結合長相互作用を使用する際は、`"Contact"`を指定します。
- `potential`: 文字列型
  - [`"Gaussian"`](GaussianPotential.md): ガウシアン型のポテンシャルを用います。
  - [`"GoContact"`](GoContactPotential.md): L-J(10,12)型のポテンシャルを用います。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `margin`: 浮動小数点型
  - カットオフ距離に対する倍率としてマージンの長さを指定します。
    定義されているポテンシャルの中で最も長いカットオフ距離に対する長さで統一されます。
