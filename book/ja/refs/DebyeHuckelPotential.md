# DebyeHückel

デバイ・ヒュッケルの式に基づいた、遮蔽効果を考慮した静電相互作用です。
溶媒を陽に扱わないシミュレーションで便利です。

{% math %}
U(r_{ij}) = \frac{q_i q_j}{4\pi\epsilon_0\epsilon_k r_{ij}} \exp(-r_{ij}/\lambda_D)
{% endmath %}

{% math %}
\lambda_D = \sqrt{\frac{\epsilon_0\epsilon_k}{2\beta N_A e_c^2 I}}
{% endmath %}

水の比誘電率の計算には、以下の論文と同じ式を使用しています。

- Sambriski, E. J. et al., (2009) Biophys. J.
- Hinckley, D. M. et al., (2013) JCP.
- Freeman, G. S., et al., (2014) JCP.

## 例

```toml
[[forcefields.global]]
# ここにはPairInteractionなどによって要求されるパラメータが入ります。
# ...
cutoff     = 5.5
parameters = [
    {index = 0, charge = 1.0},
    # 必要に応じて続きます...
]
```

## 入力

ここでの設定の他に、[`System`](System.md)の`attribute`として`temperature`と`ionic_strength`を指定する必要があります。

- `cutoff`: 浮動小数点数型
  - カットオフ距離を決めます。
  - デバイ長に対しての相対値です。書かなければ、デフォルトの値になります。
- `index`: 整数型
  - パラメータが何番目の粒子のものかを指定します。番号は0から始まります。
- `charge`: 浮動小数点数型
  - 電荷です。
