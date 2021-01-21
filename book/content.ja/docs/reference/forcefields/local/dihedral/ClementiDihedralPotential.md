+++
title = "ClementiDihedral"
weight = 400000
+++

# ClementiDihedral

ClementiらによるOff-lattice Goモデルで使用されたポテンシャルです。

{{<katex display>}}
U(v) = k_1(1-\cos(v-v_0)) + k_3(1-\cos(3(v-v_0)))
{{</katex>}}

以下の論文で提案されました。

- C. Clementi, H. Nymeyer, J. Onuchic, (2000) JMB

## 例

```toml
[[forcefields.local]]
interaction = "DihedralAngle"
potential   = "ClementiDihedral"
topology    = "none"
parameters  = [
    {indices = [0,1,2,3], v0 = -2.20, k1 = 1.0, k3 = 0.5},
    # ...
]
```

## 入力

- `v0`: 浮動小数点数型
  - ポテンシャルの最安定点を決めます。
- `k1`: 浮動小数点数型
  - このポテンシャルの強さを指定します。
- `k3`: 浮動小数点数型
  - このポテンシャルの強さを指定します。
- `indices`: 整数の配列型（長さ: 4）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。
