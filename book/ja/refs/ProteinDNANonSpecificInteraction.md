# ProteinDNANonSpecificInteraction

`ProteinDNANonSpecificInteraction`相互作用は、タンパク質-DNA間の水素結合力を模した方向依存的なコンタクトポテンシャルです。

以下の論文で提案されたポテンシャルです。

- T.Niina, G.B.Brandani, C.Tan and S.Takada, (2017) PLoS Comput Biol.

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
{index = 1000, offset = 100, kind = "Protein", PN =  999, PC = 1001, k = 1.2, r0 = 5.0, theta0 = 1.57, phi0 = 1.73},
{index = 1023, offset = 100, kind = "Protein", PN = 1022, PC = 1024, k = 1.2, r0 = 6.0, theta0 = 1.57, phi0 = 1.73},
# ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の種類を設定します。 `PDNS`です。
- `potential`: 文字列型
  - ポテンシャルの種類を設定します。`PDNS`です。
- `sigma`: 浮動小数点型
  - 距離項の届く範囲を決めるパラメータです。
- `delta`: 浮動小数点型
  - 角度項の届く範囲を決めるパラメータです。
- `parameters`: テーブルの配列型
  - `kind="DNA"`の場合、対応する3'方向の糖のインデックスが必要です。
  - `kind="Protein"`の場合、N, C末端方向の隣の粒子と、コンタクトの強さと方向のパラメータが必要です。
  - `index`: 整数型
    - どの粒子に適用するかを指定します。最初の粒子は0番目です。
  - `offset`: 整数型(optional)
    - 各粒子の粒子番号に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
