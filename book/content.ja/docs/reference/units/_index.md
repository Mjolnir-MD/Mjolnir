+++
title = "Units"
weight = 200
+++

# `[units]`
単位系を設定します。

[`[[forcefields]]`]({{<relref "/docs/reference/forcefields">}})のエネルギーのパラメータや[`[[systems]]`]({{<relref "/docs/reference/system">}})での粒子の位置など、ファイルから入力される数値は、この設定に従って使われます。

## 例

```toml
[units]
length = "angstrom"
energy = "kcal/mol"
```

## `units`

以下のフィールドがあります。

- `energy`: 文字列型
  - `"kcal/mol"`
  - `"kJ/mol"`
- `length`: 文字列型
  - `"nm"`
  - `"angstrom"`

質量の単位はAMU、電荷の単位は電気素量で固定されています。
また、角度の単位は常にラジアンです。
現在、これらは変更することはできません。

{{<hint warning>}}
時間の単位系はここで決めた単位系（と質量の単位）から計算されます。
例としては、`"kcal/mol"`と`"angstrom"`を選んだ場合、1単位時間は約49fsになります。
{{</hint>}}
