+++
title = "Contact"
weight = 40000
bookCollapseSection = true
+++

# ContactInteraction

`Contact`相互作用は、粒子間距離にかかる相互作用です。

`BondLength`相互作用と基本的に同じですが、シミュレーション中に力がかからない距離まで粒子が移動する可能性を考慮し、力がかからなくなったら力の計算をスキップするようになっています。

以下のポテンシャルが利用可能です。

- [Gaussian]({{<relref "GaussianPotential.md">}})
- [GoContact]({{<relref "GoContactPotential.md">}})
- [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
- [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})


## 例

```toml
[[forcefields.local]]
interaction = "Contact"
potential   = "GoContact"
topology    = "bond"
margin      = 0.5 # relative length to longest cutoff
parameters  = [
    # required parameters depend on a potential...
    {indices = [0, 1], ... },
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"Contact"`を指定します。
- `potential`: 文字列型
  - [Gaussian]({{<relref "GaussianPotential.md">}})
  - [GoContact]({{<relref "GoContactPotential.md">}})
  - [AttractiveGoContact]({{<relref "GoContactPotential.md">}})
  - [RepulsiveGoContact]({{<relref "GoContactPotential.md">}})
- `topology`: 文字列型
  - [Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})に設定する名前を指定します。
- `margin`: 浮動小数点型
  - カットオフ距離に対する倍率としてマージンの長さを指定します。
    定義されているポテンシャルの中で最も長いカットオフ距離に対する長さで統一されます。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型（長さ: 2）
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。
