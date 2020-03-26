# DihedralAngleInteraction

`DihedralAngle`相互作用は、二面角にかかる相互作用です。

粒子`i`番目, `j`番目, `k`番目, `l`番目の粒子について、
`i`, `j`, `k`がなす面と`j`, `k`, `l`のなす面のなす角に対してかかります。

- [Cosine](CosinePotential.md)
- [Gaussian](GaussianPotential.md)
- [ClementiDihedral](ClementiDihedralPotential.md)
- [FlexibleLocalDihedral](FlexibleLocalDihedral.md)

## 例

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Cosine"
topology    = "none"
parameters  = [
    {indices = [0, 1, 2, 3], ... }, # ポテンシャルによって要求されるパラメータは変化します。
    # 必要に応じて続きます...
]
```

## 入力

`DihedralAngle`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 二面角相互作用を使用する際は、`"DihedralAngle"`を指定します。
- `potential`: 文字列型
  - [`"ClementiDihedral"`](ClementiDihedralPotential.md): Clementi-Go力場に使用するポテンシャルです。
  - [`"Gaussian"`](GaussianPotential.md): ガウシアン型のポテンシャルを用います。
  - [`"FlexibleLocalDihedral"`](FlexibleLocalDihedral.md): AICG2+力場に使用するポテンシャルです。
  - [`"Cosine"`](CosinePotential.md): コサイン関数を使った汎用的なポテンシャルです。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `parameters`: テーブルの配列型
  - `indices`: 整数の配列型
    - どの粒子の間の距離に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - 他のパラメータは、ポテンシャルによって異なります。

## 組み合わせ

いくつかの力場では、同じ二面角に複数のポテンシャルが適用されていることがあります。
そのような場合に計算を高速化するため、以下の組み合わせがサポートされています。

それぞれのポテンシャルのパラメータは、`PotentialName = {...}`のようなインライン
テーブルとして与えられます。

{% hint style='info' %}
全ての組み合わせがサポートされているわけではありません。
{% endhint %}

### `potential = "Gaussian+FlexibleLocalDihedral"`

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Gaussian+FlexibleLocalDihedral"
topology    = "none"
env.ALA-ALA = [2.2056, 0.2183, -0.0795, 0.0451, -0.3169, 0.0165, -0.1375]
parameters  = [
{indices = [0,1,2,3], Gaussian = {v0=-2.2, k=-0.43,sigma=0.15}, FlexibleLocalDihedral = {k=1.0, coef="ALA-ALA"}},
# ...
]
```

### `potential = "Gaussian+Cosine"`

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "Gaussian+Cosine"
topology    = "none"
parameters  = [
{indices = [0,1,2,3], Gaussian = {v0=-1.57, k=-1.0,sigma=0.15}, Cosine = {k=1.0, n=1, v0 = 1.57}},
# ...
]
```

