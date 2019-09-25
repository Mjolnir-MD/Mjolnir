# DirectionalContactInteraction

`DirectionalContact`相互作用は、結合が形成される粒子間の距離とその結合方向に対してかかる相互作用です。
ポテンシャルは角度に関するポテンシャルと距離に関するポテンシャルの積として表現されます。

粒子`i`番目,`j`番目,`k`番目,`l`番目の粒子について、`i`,`j`,`k`がなす角に`angle1`で指定したポテンシャルが、`j`,`k`,`l`がなす角に`angle2`で指定したポテンシャルが、`j`,`k`間の距離に`contact`で指定したポテンシャルがかかります。

{% math %}
U(r, \theta_1, \theta_2) = U_contact (r) U_{angle_1} (\theta_1) U_{angle_2} (\theta_2)
{% math %}


角度に関する利用可能なポテンシャルには以下のようなものがあります。

- [Cosine](CosinePotential.md)

距離に関する利用可能なポテンシャルには以下のようなものがあります。

- [Gasussian](Gaussian.md)
- [GoContact](GoContact.md)

## 例

```toml
[[forcefields.local]]
interaction        = "DirectionalContact"
potentials.angle1  = "Cosine"
potentials.angle2  = "Cosine"
potentails.contact = "GoContact"
topology           = "none"
margin             = 0.5 # relative length to longest cutoff
parameters         = [
    # 各ポテンシャルそれぞれについて、要求されるパラメータを設定します。
    {indices = [0, 1, 2, 3], ... }, angle1_params = { v0 = ... }, angle2_params = { v0 = ... }, contact_params = { v0 = ... }},
    # 必要に応じて続きます...
]
```

## 入力

`DirectionalContact`相互作用を用いる場合、`[[forcefields.local]]`テーブルには、以下の値を設定します。

- `interaction`: 文字列型
  - 結合長相互作用を使用する際は、`"DirectionalContact"`を指定します。
- `potentials`
  - `angle1`,`angle2`: 文字列型
    - ["Cosine"](Cosine.md): コサイン関数を使った汎用的なポテンシャルです。
  - `contact`: 文字列型
    - ["Gaussian"](Gaussian.md): ガウシアン型のポテンシャルを用います。
    - ["GoContact"](GoContact.md): L-J(10,12)型のポテンシャルを用います。
- `topology`: 文字列型
  - [`"Topology"`](Topology.md)に設定する名前を指定します。
- `margin`: 浮動小数点型
  - `contact`で指定したポテンシャルのカットオフ距離に対する倍率としてマージンの長さを指定します。
    `contact`で定義されているポテンシャルの中で最も長いカットオフ距離に対する長さで統一されます。
