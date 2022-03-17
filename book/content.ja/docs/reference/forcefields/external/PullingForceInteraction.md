+++
title = "PullingForce"
weight = 50000
+++

# PullingForce

指定した粒子に一定の力をかけ続けます。

## 例

```toml
[[forcefields.external]]
interaction = "PullingForce"
parameters  = [
    {index = 0, force = [1.0, 0.0, 0.0]},
    {index = 0, force = 1.0, direction = [1.0, 1.0, 1.0]},
    # ...
]
```

## 入力

- `interaction`: 文字列型
  - 相互作用の名前です。ここでは、`"PullingForce"`です。
- `parameters`: テーブル型
  - `index`: 整数型
    - 力をかける粒子のインデックスです。インデックスは0始まりです。
  - `force`: 浮動小数点数型または浮動小数点数の配列型
    - 浮動小数点数の配列型が与えられた場合、力そのものを表します。
    - 浮動小数点数型が与えられた場合、力の強さを表します。その場合、`direction`を与える必要があります。
    - 単位は[`[units]`]({{<relref "../../units/_index.md" >}})で指定した単位系に従います（例えば、 `kcal/mol/Å`）.
       - `1 kcal/mol/Å ~ 69.5 pN`.
  - `direction`: 浮動小数点数の配列型
    - 力の向きです。強さは`force`で与えられるので、読み込まれたのち規格化されます。
