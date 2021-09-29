+++
title = "HybridForceField"
weight = 6000
+++

# HybridForceField

HybridForceField は二つの異なる力場の線形結合です。

{{< katex display >}}
V = \lambda V_1 + (1 - \lambda) V_2
{{< /katex >}}

トポロジーは異なっていても問題ありません。

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
