# DebyeHückel

デバイ・ヒュッケルの式に基づいた、遮蔽効果を考慮した静電相互作用です。
溶媒を陽に扱わないシミュレーションで便利です。

以下の論文で導入された力場を使用しています。

- Hinckley, D. M. et al., (2013) JCP.

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
