# GlobalStoichiometricInteraction

`GlobalStoichiometricInteraction`相互作用は、結合を形成する粒子同士が化学量論を保存するようにかかる相互作用です。

粒子`A`と`B`について、$A:B=n_a:n_b, (n_a < coefA, n_b < coefB)$の結合を許容するようなポテンシャルがかかります。

{% math %}
U(r_{ab}) = - \sum_{(a,b)} \epsilon \sigma(2 (coefB + \frac{1}{2} - \sum_{b'} u(r_{ab'})) \
                                    \sigma(2 (coefA + \frac{1}{2} - \sum_{a'} u(r_{a'b})))) u(r_{ab})
{% endmath %}

## 例

```toml
[[forcefields.global]]
interaction     = 'Stoichiometric'
potential       = 'Stoichiometric'
epsilon = 5.0
coefA = 3
coefB = 2
v0    = 4.0
range = 2.0
paramters = [
    {index = 0, kind = "A"},
    {index = 1, kind = "B"},
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - `"Stoichiometric"`を指定します。
- `potential`: 文字列型
  - 質点間の距離に適用するポテンシャルを指定します。
- `parameters`: テーブルの配列型
  - `index`: 整数型
    - 適用する粒子の番号を指定します。
  - `offset`: 整数型(optional)
    - `indices`に加算する値です。省略可能です。グループ内番号を使いたい場合に便利です。
  - `kind`: 文字列型
    - 粒子が`A`であるか`B`であるかを指定します。それ以外の文字列を指定した場合はエラーになります。
