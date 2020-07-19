+++
title = "AFMFitting"
weight = 40000
+++

# AFMFittingInteraction

AFM画像への構造のフレキシブルフィッティングを行います。

以下の論文で提案されました。

- T. Niina et al., JCTC (2020)

{{<katex display>}}
\begin{aligned}
U(\mathbf{r}) &= k(1 - \mathrm{c.c.}(\mathbf{r})) \\
\mathrm{c.c.} &= \frac{\sum_{p\in\mathrm{pixels}} H_p^{\mathrm{(exp)}} H_p^{\mathrm{(sim)}}(\mathbf{r})}
                     {\sqrt{\sum_{p\in\mathrm{pixels}} \left(H_p^{\mathrm{(exp)}}\right)^2}
                      \sqrt{\sum_{p\in\mathrm{pixels}} \left(H_p^{\mathrm{(sim)}}(\mathbf{r})\right)^2}} \\
H_p^{\mathrm{(sim)}}(\mathbf{r}) &= \gamma\log\left(1 + \sum_i^N \exp\left(\frac{-(x_i - x_p)^2 - (y_i - y_p)^2 }{2\sigma^2}\right)\exp\left(\frac{z_i + r_i}{\gamma}\right)\right)
\end{aligned}
{{</katex>}}

## 例

```toml
[[forcefields.external]]
interaction = "AFMFlexibleFitting"
k           = 100.0
gamma       =   1.0
pixel_x     =  10.0
pixel_y     =  10.0
length_x    =   5
length_y    =   5
sigma_x     =   2.0
sigma_y     =   2.0
z0          = 0.0
cutoff      = 5.0
margin      = 0.5
image       = [
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, 0.0,
    0.0, 0.0, 0.5, 1.0, 0.5,
    0.0, 0.0, 1.0, 2.0, 1.0,
    0.0, 0.0, 0.5, 1.0, 0.5,
]
parameters  = [
{index = 0, radius = 1.0},
{index = 1, radius = 2.0},
{index = 4, radius = 3.0},
{index = 5, radius = 4.0},
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。ここでは`"AFMFlexibleFitting"`です。
- `k`: 浮動小数点数型
  - ポテンシャルの強さを決めます。
- `gamma`: 浮動小数点数型
  - 関数{{<katex>}}H{{</katex>}}のパラメータで、ソフトマックスの正確さを決めます。
- `pixel_x, pixel_y`: 浮動小数点数型
  - それぞれの方向でのピクセルの大きさです。
- `length_x, length_y`: 整数型
  - それぞれの方向でのピクセルの数です。
- `sigma_x, sigma_y`: 浮動小数点数型
  - 関数{{<katex>}} H {{</katex>}} のパラメータで、ピクセルの影響が届く範囲を決めます。
- `z0`: 浮動小数点数型
  - 数値誤差を抑えるための内部オフセットです。通常0で構いません。
- `cutoff`: 浮動小数点数型
  - ピクセルの影響をカットオフする長さです。{{<katex>}} \sigma {{</katex>}}との相対です。
- `margin`: 浮動小数点数型
  - 内部の近接リストで用いるマージンです。カットオフに対する相対です。
- `image`: 浮動小数点数の配列型
  - 参照AFM画像です。各ピクセルの高さを設定します。
  - 1つめの要素がピクセル (0, 0)を、2つめが (1, 0)、... (Lx, 0)、(0, 1) pixel, ... と続きます。
  - (0, 0)ピクセルは`(0.0, 0.0)`から`(pixel_x, pixel_y)`までの範囲を占める四角形です。
  - (n, m)ピクセルは`(n*pixel_x, m*pixel_y)`から`((n+1)*pixel_x, (m+1)*pixel_y)`までの範囲を占める四角形です。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 粒子の番号です。最初の粒子は0番目です。
  - `radius`: 浮動小数点数型
    - 粒子のAFM観測下での半径です。
