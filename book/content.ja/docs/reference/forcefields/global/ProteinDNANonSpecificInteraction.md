+++
title = "PDNS"
weight = 30000
+++

# ProteinDNANonSpecificInteraction

`ProteinDNANonSpecificInteraction`は粗視化モデルにおいてDNAとタンパク質の間の水素結合をモデル化したものです。

以下の論文で導入されました。

- T.Niina‡, G.B.Brandani‡, C.Tan and S.Takada, (2017) PLoS Comput Biol. (‡: co-1st)

## 例

```toml
[[forcefields.global]]
interaction = "PDNS"
potential   = "PDNS"
spatial_partition.type   = "VerletList"
spatial_partition.margin = 0.5
sigma = 1.0
delta = 0.17453
parameters  = [
{index =    2, kind = "DNA", S3 = 1},
{index =    5, kind = "DNA", S3 = 4},
# ...
{index = 1000, offset = 100, kind = "Protein", PN =  999, PC = 1001, k = -1.2, r0 = 5.0, theta0 = 1.57, phi0 = 1.73},
{index = 1023, offset = 100, kind = "Protein", PN = 1022, PC = 1024, k = -1.2, r0 = 6.0, theta0 = 1.57, phi0 = 1.73},
# ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。`"PDNS"`です。
- `potential`: 文字列型
  - ポテンシャルの名前です。`"PDNS"`です。
- `sigma`: 浮動小数点数型
  - 動径方向での安定化領域の幅です。
- `delta`: 浮動小数点数型
  - 角度方向での安定化領域の幅です。単位はラジアンです。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 粒子の番号です。
    - `DNA`の場合、リン酸残基の番号です。
  - `offset`: 整数型（省略可能）
    - 粒子のオフセットです。
  - `kind`: 文字列型
    - `DNA`粒子と`Protein`粒子があります。
  - `S3`: 整数型
    - `DNA`粒子で要求されます。3'末端方向にある隣のヌクレオチドの糖粒子の番号です。
  - `PN, PC`: 整数型
    - `Protein`粒子で要求されます。それぞれN-, C-末端方向の隣の粒子の番号です。
  - `k`: 浮動小数点数型
    - ポテンシャルの強さです。
  - `r0`, `theta0`, `phi0`: 浮動小数点数型
    - 再安定構造です。
