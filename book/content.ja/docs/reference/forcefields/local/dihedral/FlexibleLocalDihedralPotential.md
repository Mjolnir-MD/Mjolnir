+++
title = "FlexibleLocal"
weight = 300000
+++

# FlexibleLocalDihedral

粗視化タンパク質モデルで、フレキシブル領域の角度分布を再現するためのポテンシャルです。
角度分布を再現するため、フーリエ級数展開の形で定義されています。

{{<katex display>}}
U(\phi) = C + \sum_{n=1}^{3}\left( k_n^{\sin} \sin(n\phi) + k_n^{\cos} \cos(n\phi)\right)
{{</katex>}}

AICG2+粗視化タンパク質力場の一部として使われることが多いです。

以下の論文で開発されました。

- T. Terakawa and S. Takada, (2011) Biophys J 

## 例

```toml
[[forcefields.local]]
parameters = [
    {indices = [0,1,2,3], k = 1.0, coef = [2.2356,  0.4119, -0.1283,  0.0229, -0.2708, -0.0085, -0.0641]},
    # ...
]
```

## 入力

- `k`: 浮動小数点数型
  - このポテンシャルの強さを指定します。
- `coef`: 浮動小数点数の配列型（長さ: 7）
  - 前から順に、{{<katex>}} C, k_1^{\sin}, k_1^{\cos}, k_2^{\sin}, k_2^{\cos}, k_3^{\sin}, k_3^{\cos} {{</katex>}}です。
- `indices`: 整数の配列型（長さ: 4）
  - どの粒子の間に適用するかを指定します。最初の粒子は0番めです。
- `offset`: 整数型（省略可能）
  - インデックスに加算する値です。省略可能です。

