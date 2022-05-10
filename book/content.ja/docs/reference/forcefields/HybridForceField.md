+++
title = "HybridForceField"
weight = 6000
+++

# HybridForceField

HybridForceField は二つの異なる力場の線形結合です。

{{< katex display >}}
V = \lambda V_1 + (1 - \lambda) V_2
{{< /katex >}}

{{< katex >}}\lambda = 1{{< /katex >}}のとき、{{< katex >}}V = V_1{{< /katex >}}になることに気を付けてください。

トポロジーは異なっていても問題ありません。異なっている場合、dcdなどのchain数には一つ目のforcefieldが入ります。

## 例

まず、二つの`forcefields`を定義します。

その後、`[simulator]`テーブルで`forcefields.type = "Hybrid"`と、`forcefields.lambda`の値を定義してください。
一番目の`[[forcefields]]`テーブルが上の式での{{<katex>}} V_1 {{</katex>}}に、二番目が{{<katex>}} V_2 {{</katex>}}になります。

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 1000000
save_step     =   1_000
seed          = 2859805901
integrator.type = "BAOABLangevin"
integrator.gammas = [
    # ...
]
forcefields.type = "Hybrid"
forcefields.lambda = 0.5
# ...

[[forcefields]] # V1
[[forcefields.local]]
# ...

[[forcefields]] # V2
[[forcefields.local]]
# ...
```

## λの値を動的に変化させる

{{< katex >}}\lambda{{< /katex >}}を仮想的な粒子の座標と解釈し、それにかかる力を考え、時間積分することによってサンプリングを行うことが可能です。
`forcefields.lambda`には文字列で`"dynamic"`と入力し、`[[systems]]`で`dynamic_variables`に`lambda`を定義してください。
この`lambda`の値、速度、力は`.ene`ファイルに保存されます。

また、その際に{{< katex >}}k(\lambda - \lambda_0)^2{{< /katex >}}型のポテンシャルを追加でかけることが可能です。

```toml
[simulator]
type          = "MolecularDynamics"
boundary_type = "Unlimited"
precision     = "double"
delta_t       = 0.1
total_step    = 1000000
save_step     =   1_000
seed          = 2859805901
integrator.type = "BAOABLangevin"
integrator.gammas = [
    # ...
]
forcefields.type = "Hybrid"
forcefields.lambda = "dynamic" # make lambda dynamic
forcefields.v0_lambda = 1.0    # k(v - v0)^2
forcefields.k_lambda  = 100.0  # k(v - v0)^2

# ...

[[systems]]
attributes.temperature = 300.0
dynamic_variables.lambda = {x = 1.0, m = 1000.0, gamma = 1e-5, boundary = "Repulsive", lower = 0.0, upper = 1.0}

[[forcefields]] # V1
[[forcefields.local]]
# ...

[[forcefields]] # V2
[[forcefields.local]]
# ...
```
