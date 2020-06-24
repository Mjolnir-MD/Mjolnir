# MultipleBasinForceField

MultipleBasin モデルによって、複数の異なる力場をスムースにつなげることができます。

粗視化モデルでタンパク質の構造変化を実現するための、力場の強制変更よりも自然な方法として、以下の文献で提案されました。

- Kei-ichi Okazaki, Nobuyasu Koga, Shoji Takada, Jose N. Onuchic, and Peter G. Wolynes PNAS (2006)

Mjolnirはこれを一般化した形で実装しており、任意の力場を利用できます。

このモデルでは、以下の形で定義された行列の最小固有値として定義される $$ V_{MB} $$ がかかります。

{% math %}
\begin{pmatrix}
V_1 + \Delta V_1 & \Delta \\
\Delta & V_1 + \Delta V_1
\end{pmatrix}
\begin{pmatrix}
c_1 \\ c_2
\end{pmatrix}
= V_{MB}
\begin{pmatrix}
c_1 \\ c_2
\end{pmatrix}
{% endmath %}

固有ベクトルの成分は、それぞれの状態の重みと解釈できます。反応座標として、デフォルトで以下の値が出力されます。

{% math %}
\chi = \log\left(\frac{c_2}{c_1}\right)
{% endmath %}

$$ 0 < \Delta $$ のとき、この値は非数になることに注意してください。
これを避けるため、Mjolnirでは $$\Delta$$ を常に負の値として読み込みます。
$$ \Delta $$ は $$ V_{MB} $$の値には二乗の形でしか効かないため、その絶対値にのみ意味があります。

3つの異なる力場をつなげる際には、以下の形になります。

{% math %}
\begin{pmatrix}
V_1 + \Delta V_1 & \Delta_{12} & \Delta_{13} \\
\Delta_{21} & V_2 + \Delta V_2 & \Delta_{23} \\
\Delta_{31} & \Delta_{32} & V_3 + \Delta V_3 \\
\end{pmatrix}
\begin{pmatrix}
c_1 \\ c_2 \\ c_3
\end{pmatrix}
= V_{MB}
\begin{pmatrix}
c_1 \\ c_2 \\ c_3
\end{pmatrix}
{% endmath %}

ここで、$$\Delta_{ij} = \Delta_{ji}$$ です。

現在、Mjolnirは力場の数が2つまたは3つの場合のみをサポートしています。
さらに、異なるトポロジーを持つ力場同士のMultipleBasinはサポートしていません。

## 例

まず、繋ぎ合わせるための複数種類の力場を設定し、名前をつけます。
この時、`"common"`という名前の力場はMultipleBasinに含まれず常に適用されるようになることに注意してください。

その後、`[simulator]`でMultipleBasinを利用することを宣言し、力場の名前と、カップリングパラメータを指定します。

独立に構造変化を起こす複数のタンパク質をサポートする場合などのために、複数のMultipleBasinを並行して定義できるようになっています。
`simulator.forcefields.units`で必要な数のMultipleBasinを定義してください。

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
forcefields.type = "MultipleBasin"
forcefields.units = [
{basins = ["open", "close"], dVs = [0.0, -12.0], delta = 150.0},
]
# ...

[[forcefields]]
name = "open"
[[forcefields.local]]
# ...

[[forcefields]]
name = "close"
[[forcefields.local]]
# ...

[[forcefields]]
name = "common"
[[forcefields.local]]
# ...
```

## AICG2+力場で用いる時の注意

AICG2+力場を用いる際は、結合長とコンタクトの扱いに注意が必要です。
遷移状態のエネルギーを上げすぎないように、自然長が大きく違う場合、以下の形でポテンシャルの係数を補正します。

{% math %}
K'_b = K_b \times \mathrm{min}\left(1, \frac{E_{max}}{K_b(b_i^{(1)} - b_i^{(2)})^2}\right)
{% endmath %}

また、コンタクトは引力項と斥力項に分割され、最も短い参照距離を持つコンタクトの斥力項を全ての力場で共有します。
このときコンタクトがもともと存在しない力場においては、斥力項のみが与えられ、代わりに排除体積効果は無視されます。
引力項では、その力場がもともと持っていたコンタクトのパラメータを使用します。
