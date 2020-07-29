+++
title = "Dummy"
weight = 60000
+++

# DummyInteraction

`Dummy`相互作用は、何も計算しない相互作用です。利用可能なポテンシャルはありません。

これは、力場に必要な相互作用だけでは表現できないようなトポロジーを表現するために
用いることができます。つまり、これを使うと、実際には直接の相互作用がないような
粒子のペアに対して結合を加え、それを基準にnonbondedな相互作用に影響を与えること
ができます。

トポロジーを指定しなかった（`"none"`にした）場合、この相互作用には何の意味も
なくなります。

トポロジーの効果については、[Topology]({{<relref "/docs/reference/forcefields/Topology.md">}})を参照して下さい。

## 例

```toml
[[forcefields.local]]
interaction = "Dummy"
# No potential field required.
topology    = "bond"
parameters  = [
    {indices = [0, 1], offset = 100}, # No other parameters.
    # ...
]
```

## 入力

- `interaction`: String
  - 相互作用の名前です。`"Dummy"`を指定します。
- `topology`: 文字列型
  - [`"Topology"`]({{<relref "/docs/reference/forcefields/Topology.md">}})に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型（長さ: 2）
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(省略可能)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - ポテンシャルがないので、他のパラメータはありません。
